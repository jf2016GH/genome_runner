import cherrypy, os, cgi, tempfile, sys, itertools
from mako.lookup import TemplateLookup
from mako.exceptions import RichTraceback as MakoTraceback
from path import PathNode
import path as grsnp_path
from operator import itemgetter
from path import base_name as basename
import logging
from logging import FileHandler,StreamHandler
import json
import pdb
import grsnp.worker_gr
from time import gmtime, strftime
import simplejson
import string
import random
import traceback
import dbcreator_ucsc as uscsreader
import argparse
import shutil
import subprocess
from celery import Celery
import server_utils as utils


os.environ['GR_COMPATIBILITY_MODE'] = 'y'

root_dir = os.path.dirname(os.path.realpath(__file__))
lookup = TemplateLookup(directories=[os.path.join(root_dir,"frontend/templates")])


sett = {}
DEBUG_MODE = True
logger = logging.getLogger('genomerunner.server')





# Each function in this class is a web page 
class WebUI(object):
	def __init__(self):		
		# go through each database directory and create custom_data if it does not exist.
		for db_ver,db_dir in sett["data_dir"].items():
			# create all directories in the custom_data dir if they do not already exist
			for org in self.get_org(db_ver):
				custom_dir = os.path.join(os.path.split(db_dir)[0],"custom_data")
				logger.info("Processing genomic features for {}".format(org))
				if not os.path.exists(custom_dir): os.mkdir(custom_dir)
				cust_sub_dir = ["backgrounds","gfs","fois","rsid_conversion"]
				for c in cust_sub_dir:
					tmp = os.path.join(custom_dir,c)
					if not os.path.exists(tmp): os.mkdir(tmp)
					c_dir = os.path.join(custom_dir,c,org)
					if not os.path.exists(c_dir): os.mkdir(c_dir)
				# Read the genomic feature files and generate html files
				paths = PathNode()
				paths.name = "Root"
				paths.organisms = self.get_org(db_ver) 
				paths.traverse(os.path.join(db_dir,org))
				grsnp_path.write_treeview_json(os.path.join(db_dir,org))
		self._index_html = {}

	@cherrypy.expose
	def index(self,organism=None,db_version=None):
		global results_dir, uploads_dir, sett
		if not organism: organism = sett["default_organism"]

		if DEBUG_MODE or not organism in self._index_html:		
			tmpl = lookup.get_template("master.mako")
			paths = PathNode()
			paths.name = "Root"
			[html_dbversion, db_version] =  grsnp_path.get_database_versions_html(sett["data_dir"],db_version)
			paths.organisms = self.get_org(db_version)
			# Check if the organism actually exists in the current database.
			# If it is not, select the default organism
			if not organism in paths.organisms:
				organism = sett["default_organism"]
			custom_dir = os.path.join(os.path.split(sett["data_dir"][db_version])[0],"custom_data")
			# Use mako to render index.html
			body = lookup.get_template("index.mako").render(paths=paths,default_background=paths.get_backgrounds_combo(organism,custom_dir),
									custom_gfs=paths.get_custom_gfs(organism,custom_dir),demo_snps=paths.get_custom_fois(organism,custom_dir),
									data_dir=os.path.join(sett["data_dir"][db_version],organism),default_organism=organism,
									database_versions=html_dbversion,pct_scores=paths.get_scores(os.path.split(sett["data_dir"][db_version])[0]))
			script = lookup.get_template("index.js").render(default_organism=organism)
			self._index_html[organism] = tmpl.render(body=body,script=script)

		return self._index_html[organism]

	def get_org(self,db_version):
		organisms = []
		files = os.listdir(sett["data_dir"][db_version])
		for f in files:
			if f.find(".") == -1:
				organisms.append(f)
		return organisms	

	@cherrypy.expose
	def query(self, strand="",run_annotation=False,db_version=None,jstree_gfs="",stat_test=None,num_mc = None,**kwargs):
		global results_dir, uploads_dir, sett
		# Assign a random id
		id = ''.join(random.choice(string.lowercase+string.digits) for _ in range(32))
		while (os.path.exists(os.path.join(uploads_dir,id))):
			id = ''.join(random.choice(string.lowercase+string.digits) for _ in range(32))
		res_dir = os.path.join(results_dir,str(id))
		upload_dir = os.path.join(uploads_dir,str(id))
		os.mkdir(upload_dir)
		os.mkdir(os.path.join(upload_dir,"fois"))
		os.mkdir(os.path.join(upload_dir,"gfs"))
		res_dir = os.path.join(results_dir,str(id))
		os.mkdir(res_dir)

		fois = os.path.join(upload_dir,".fois") # contains a list of the paths to fois to run through the analysis
		gfs = os.path.join(upload_dir,".gfs") # contains a list of the paths to the gfs to run the fois against
		list_gfs = []
		data_dir = os.path.split(sett["data_dir"][db_version])[0]
		cherrypy.response.timeout = 3600
		# fill in missing keys with "" for input controls that are blank
		required_keys = ["bedfile:fois", "bedcustom:fois", "bedtext:fois", "bedfiles:gfs", "bedtext:gfs",
						 "bedcustom:gfs", "bedfile:background", "bedtext:background","bedpath:background"]
		for r in required_keys:
			if r not in kwargs.keys():
				kwargs[r] = ""


		# get the FOI files, group box selection, and text that are uploaded by the user
		list_fois_paths = []  # contains paths to all uploaded foi data
		fois_outdir = os.path.join(upload_dir,"fois")
		list_fois_paths = list_fois_paths + utils.retrieve_files(kwargs["bedfile:fois"],fois_outdir,id)
		list_fois_paths = list_fois_paths + utils.retrieve_group(kwargs["bedcustom:fois"])
		# only use the textbox if nothing been uploaded
		if len(list_fois_paths) == 0:
			list_fois_paths = list_fois_paths + utils.retrieve_text(kwargs["bedtext:fois"],os.path.join(fois_outdir,"custom_foi.bed"),id)
		# remove fois with the same base_name
		purged_fois_paths = utils.purge_duplicate_features(list_fois_paths)
		# write foi paths to .fois file
		with open(fois, 'a') as out_fois:
			for f in purged_fois_paths:
				out_fois.write(f + "\n")
		if len(purged_fois_paths) == 0:
			return "Feature of Interest files not detected.  Please upload or choose Feature of Interest to run."

		# get the GFs
		list_gfs_paths = []
		gfs_outdir = os.path.join(upload_dir, "gfs")
		# Create annotation folder, used by the server to check if annotation is going to be run
		list_gfs_paths = list_gfs_paths + utils.retrieve_files(kwargs["bedfile:gfs"],gfs_outdir,id)
		# collect the gfs for all groups that were checked by the user
		checked_custom_gfs = { k: kwargs[k] for k in kwargs.keys() if k.startswith("bedcustom:gfs") and kwargs[k] == "on"}
		for k,v in checked_custom_gfs.iteritems():
			gp_gfs_dir = k.split(":")[-1]
			list_gfs_paths = list_gfs_paths + utils.retrieve_group(gp_gfs_dir)

		# add genomic feature tracks loaded via the JStree control
		for k in jstree_gfs.split(','):
			if k.startswith('file:'):
				g = k.split(":")[-1]
				list_gfs_paths.append(g)
		# only use the textbox if nothing been uploaded
		if len(list_gfs_paths) == 0:
			list_gfs_paths = list_gfs_paths + utils.retrieve_text(kwargs["bedtext:gfs"], os.path.join(gfs_outdir,"custom_gf.bed"), id)
		# remove gfs with the same base_name
		purged_gfs_paths = utils.purge_duplicate_features(list_gfs_paths)
		# write all gfs paths to .gfs file
		with open(gfs, "a") as out_gfs:
			for g in purged_gfs_paths:
				tmp_g = verify_score_strand(g, kwargs['pct_score'], strand, data_dir)
				out_gfs.write(tmp_g + "\n")

		if run_annotation:
			annot_outdir = os.path.join(res_dir,"annotations")
			if not os.path.exists(annot_outdir): os.mkdir(annot_outdir)


		# load the background data if uploaded.
		background_path = utils.retrieve_files(kwargs["bedfile:background"], upload_dir,id)
		if len(background_path) == 0:
			background_path = utils.retrieve_text(kwargs["bedtext:background"], os.path.join(upload_dir,"custom.background.bed"), id)
		if len(background_path) !=0: background_path = background_path[0]
		else:
			background_path = kwargs["bedpath:background"] # returns a string

		# make paths relative, needed for remote celery workers to function correctly
		for f in [fois,gfs]:
			list_foi = open(f).read().replace(uploads_dir,'/uploads').replace(results_dir,'/results').replace(data_dir,"")
			with open(f,'wb') as writer:
				writer.write(list_foi)
				background_path = background_path.replace(data_dir,'').replace(os.path.split(uploads_dir)[0],"").lstrip("/")

		if stat_test == "montecarlo":
			stat_test = "montecarlo_"+num_mc
		organism = "",False
		for k,v in kwargs.items():
			# organism to use
			if "organism:" in v:
				organism = v.split(":")[-1]
			if k.startswith("run_annot") and v == "on": run_annotation = True
		# write the enrichment settings.
		path = os.path.join(res_dir, ".settings")
		set_info = {"Jobname:": str(id),
					"Time:": strftime("%Y-%m-%d %H:%M:%S", gmtime()),
					"Background:": os.path.split(background_path)[-1],
					"Organism:": organism,
					"Database version:":db_version,
					"% Score threshold:": str(kwargs['pct_score'])+"%",
					"Strand:": strand,
					"P-value statistical test:": stat_test}

		with open(path, 'wb') as sett_files:
			for k,v in set_info.iteritems():
				sett_files.write(k+"\t"+v+"\n")

		gfs_count = 0
		if open(gfs).read() == "":
			return "ERROR: No Genomic Features selected/uploaded."
		else:
			gfs_count = len([x for x in open(gfs).read().split("\n") if x != ""])

		# copy over .foi to results folder, used for the enrichment results
		shutil.copy(fois,os.path.join(res_dir,".fois"))

		# run using celery queues.
		run_args = [fois,gfs,background_path,id,True,os.path.join(sett["data_dir"][db_version],organism,"bkg_overlaps.gr"),run_annotation,False]
		run_kwargs = { "pct_score": kwargs['pct_score'],"organism": organism,"id": id,"db_version": db_version,"stat_test": stat_test }
		if gfs_count > 3:
			print "LONG RUN STARTED"
			run_queue = 'long_runs'
		else:
			print "SHORT RUN STARTED"
			run_queue = 'short_runs'
		try:
			grsnp.worker_gr.run_hypergeom.apply_async(args=run_args, kwargs = run_kwargs, queue=run_queue, retry=False)
		except Exception, e:
			print "WORKER ERROR"
		raise cherrypy.HTTPRedirect("result?id=%s" % id)

	@cherrypy.expose
	def result(self, id):
		global results_dir, uploads_dir, sett
		path = os.path.join(results_dir, id)
		params = {}
		params["run_id"] = id
		params["detailed"] = "Results not yet available"
		params["matrix"] = "Results not yet available"		
		tmpl = lookup.get_template("master.mako")

		# Loads the progress file if it exists
		p = {"status":"","curprog":0,"progmax":0}
		progress_path = os.path.join(path,".prog")
		if os.path.exists(progress_path):
			with open(progress_path) as f:
				p = json.loads(f.read() )
		params["log"] = "###Run Settings###\n"
		sett_path = os.path.join(path,".settings")
		organism = ""
		if os.path.exists(sett_path):
			with open(sett_path) as f:
				tmp = f.read()	
				params["log"] = params["log"]+ tmp
				organism = [x.split("\t")[1] for x in tmp.split("\n") if x.split("\t")[0] == "Organism:"][0]
		params["organism"] = organism
		params["log"] = params["log"] + "\n###Run Log###\n"
		debug_path = os.path.join(path,".log")
		if os.path.exists(debug_path):
			with open(debug_path) as f:
				params["log"] = params["log"] + f.read()

		# loads results from results file		
		detailed_path = os.path.join(path,"detailed.gr")
		if os.path.exists(detailed_path):
			with open(detailed_path) as f:
				params["detailed"] = f.read()
		
		foi_names_path = os.path.join(os.path.join(results_dir, id),".fois")
		if os.path.exists(foi_names_path):
			with open(foi_names_path) as f:
				params["fois"] = [basename(x).split(".")[0] for x in f.read().split("\n") if x != ""]
		else:
			params["fois"] = ""

		params["zipfile"] = os.path.join("results",id,"GR_{}.tar.gz").format(id)
		params["run_annotation"] = True if os.path.exists(os.path.join(results_dir,id,"annotations")) else  False
		params.update(p)
		try:
			rend_template = tmpl.render(body=lookup.get_template("results.mako").render(**params),script= lookup.get_template("results.js").render(**params))
			print "LOADED TEMPLATE"
		except Exception, e:
			traceback = MakoTraceback()
			str_error = ""
			for (filename, lineno, function, line) in traceback.traceback:
				str_error +=  "File %s, line %s, in %s" % (os.path.split(filename)[-1], lineno, function)
				str_error += "\n"
				str_error += line + "\n"
				str_error += "%s: %s" % (str(traceback.error.__class__.__name__), traceback.error)
			print str_error
			rend_template = str_error
		return rend_template

	@cherrypy.expose
	def results_shiny(self, id):
		global results_dir, uploads_dir, sett
		path = os.path.join(results_dir, id)	
		params = {}	
		params['run_id'] = id
		try:
			tmp = lookup.get_template("master.mako")
			script = lookup.get_template("results_shiny.js").render(run_id=id)
			rend_template = lookup.get_template("results_shiny.mako").render(run_id=id,script=script)

		except Exception, e:
			traceback = MakoTraceback()
			str_error = ""
			for (filename, lineno, function, line) in traceback.traceback:
				str_error +=  "File %s, line %s, in %s" % (os.path.split(filename)[-1], lineno, function)
				str_error += "\n"
				str_error += line + "\n"
				str_error += "%s: %s" % (str(traceback.error.__class__.__name__), traceback.error)
			print str_error
			rend_template = str_error
		return rend_template

	@cherrypy.expose
	def gf_descriptions(self,db_version,organism):
		# Use mako to render index.html
		tmpl = lookup.get_template("master.mako")
		body = lookup.get_template("gf_descriptions.mako").render()
		script = lookup.get_template("gf_descriptions.js").render(db_version=db_version,organism=organism)
		return tmpl.render(body=body,script=script)

	@cherrypy.expose
	def get_detailed(self,run_id):
		""" loads results from detailed results file
		"""
		global results_dir, uploads_dir, sett
		detailed_path,results = os.path.join(results_dir, run_id,"detailed.txt"),{"detailed": ""}		 
		if os.path.exists(detailed_path):
			with open(detailed_path) as f:
				results["detailed"] = f.read()
		return simplejson.dumps(results)

	@cherrypy.expose
	def get_progress(self, run_id):
		# Loads the progress file if it exists
		global results_dir, uploads_dir, sett
		p = {"status":"","curprog":0,"progmax":0}
		progress_path = os.path.join(os.path.join(results_dir, run_id),".prog")
		if os.path.exists(progress_path):
			with open(progress_path) as f:
				p = json.loads(f.read())
		return simplejson.dumps(p)

	@cherrypy.expose
	def get_log(self,run_id):
		global results_dir, uploads_dir, sett
		results = {"log": ""}
		log_path = os.path.join(os.path.join(results_dir, run_id),"gr_log.txt")
		if os.path.exists(log_path):
			with open(log_path) as f:
				results["log"] = f.read()
		return simplejson.dumps(results)

	@cherrypy.expose
	def get_gfdescriptions(self,organism,db_version):
		return open(os.path.join(sett["data_dir"][db_version],organism,"gf_descriptions.txt")).read()
	@cherrypy.expose
	def get_checkboxtree(self,organism,db_version):
		return open(os.path.join(sett["data_dir"][db_version],organism,"treeview.json")).read()


	@cherrypy.expose
	def enrichment_log(self, id):
		global results_dir
		with open(os.path.join(results_dir,id+".log")) as sr:
			x = sr.read()
			return "<p>{}</p>".format(x.replace("\n","<br/>"))

	@cherrypy.expose
	def cite(self):
		tmpl = lookup.get_template("master.mako")
		return tmpl.render(body=lookup.get_template("cite.mako").render(),script="")

	@cherrypy.expose
	def news(self):
		tmpl = lookup.get_template("master.mako")
		return tmpl.render(body=lookup.get_template("news.mako").render(),script="")

	@cherrypy.expose
	def overview(self):
		tmpl = lookup.get_template("master.mako")
		return tmpl.render(body=lookup.get_template("overview.mako").render(),script="")

	@cherrypy.expose
	def demo(self):
		tmpl = lookup.get_template("master.mako")
		return tmpl.render(body=lookup.get_template("demo.mako").render(),script="")

	@cherrypy.expose
	def roadmap(self):
		tmpl = lookup.get_template("master.mako")
		return tmpl.render(body=lookup.get_template("roadmap.mako").render(),script="")

	@cherrypy.expose
	def help(self):
		tmpl = lookup.get_template("master.mako")
		return tmpl.render(body=lookup.get_template("help.mako").render(),script="")

def base_name(k):
    return os.path.basename(k).split(".")[0]

def verify_score_strand(gf_path,pct_score,strand,data_dir):
    ''' Checks if a score and/or strand filtered version of gf_path exists in the database and
    returns the appropriate path if it does.
    '''
    global results_dir, uploads_dir, sett
    gf_path = os.path.join(data_dir,gf_path.lstrip("/"))
    gf_score_strand_path = gf_path.replace('/grsnp_db/','/grsnp_db_{}_{}/'.format(pct_score,strand))
    gf_score_path = gf_path.replace('/grsnp_db/','/grsnp_db_{}/'.format(pct_score))
    gf_strand_path = gf_path.replace('/grsnp_db/','/grsnp_db_{}/'.format(strand))
    if os.path.exists(gf_score_strand_path):
        return gf_score_strand_path
    elif os.path.exists(gf_score_path):
        return gf_score_path
    elif os.path.exists(gf_strand_path):
    	return gf_strand_path
    else:
    	return gf_path

def main():
	global sett, results_dir, uploads_dir
	root_dir = os.path.dirname(os.path.realpath(__file__))
	static_dir = os.path.abspath(os.path.join(root_dir, "frontend/static"))
	media = os.path.abspath(os.path.join(".","frontend/media"))
	parser = argparse.ArgumentParser(prog="python -m grsnp.server", description="Starts the GenomeRunner SNP server. Example: python -m grsnp.server -d /home/username/db_#.##_#.##.####/ -g hg19 -p 8000", epilog="Use GenomeRunner SNP: http://localhost:8000/gr")
	parser.add_argument("--data_dir" , "-d", nargs="?",type=str, help="Set the directory containing the database. Required. Use absolute path. Example: /home/username/db_#.##_#.##.####/.", required=True)
	parser.add_argument("--run_files_dir" , "-r", nargs="?", help="Set the directory where the server should save results. Required. Use absolute path. Example: /home/username/run_files/.", required=True)
	parser.add_argument("--organism" , "-g", nargs="?", help="The UCSC code for the organism to use. Default: hg19 (human). Data for the organism must exist in the database directory. Use dbcreator to make the database, if needed.", default="hg19")
	parser.add_argument("--port","-p", nargs="?", help="Socket port to start server on. Default: 8000", default=8080) 
	parser.add_argument("--num_workers", "-w", type=int, help="The number of worker processes to start. If celery worker already exists, this is ignored. Default: 1", default=1)
	parser.add_argument("--group", "-z", type=str, help="The group to change results folder permission to", default="")	

	args = vars(parser.parse_args())
	port = args["port"]
	list_data_dir = args["data_dir"].split(",")

	if list_data_dir == "":
		print "ERROR: data_dir is a required argument"
		sys.exit()

	if args['run_files_dir'] == "":
		print "ERROR: run_files_dir is a required argument"
		sys.exit()
	data_dir = {}
	for db_dir in list_data_dir:
		if db_dir.strip()[-1] == "/":
			data_dir.update({os.path.split(db_dir[:-1])[1]:os.path.join(db_dir.strip(),"grsnp_db")})
		else:
			data_dir.update({os.path.split(db_dir)[1]:os.path.join(db_dir.strip(),"grsnp_db")})

	# setup the logging
	hdlr = logging.FileHandler(os.path.join(args['run_files_dir'], 'genomerunner_server.log'))
	hdlr_std = StreamHandler()
	formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
	hdlr.setFormatter(formatter)
	logger.addHandler(hdlr)
	logger.addHandler(hdlr_std)
	logger.setLevel(logging.INFO)

	# global settings used by GR
	sett = {"data_dir":data_dir,
				"run_files_dir": args["run_files_dir"],
				"default_organism":args["organism"],
				"group": args['group']	
				}
	
	#validate data directory
	for k,v in data_dir.items():
		if not os.path.exists(v):
			print "ERROR: {} does not exist. Please run grsnp.dbcreator or use --data_dir.".format(v)
			sys.exit()
		if not os.path.exists(os.path.join(v,sett["default_organism"])):
			print "ERROR: Database for default organism {} does not exist. Either change the default organism or add data for that organism to the database at {} using the bedfilecreator".format(sett["default_organism"],v)
			sys.exit()

	# validate run_files directory
	if not os.path.exists(sett["run_files_dir"]): os.mkdir(sett["run_files_dir"])
	results_dir = os.path.join(sett["run_files_dir"],"results")
	uploads_dir = os.path.join(sett["run_files_dir"],"uploads")	
	if not os.path.exists(results_dir):
		os.mkdir(results_dir)
	if not os.path.exists(uploads_dir):
		os.mkdir(uploads_dir)
	if port:
		cherrypy.server.max_request_body_size = 0
		cherrypy.config.update({
			"server.socket_port": int(port),
			"server.socket_host":"0.0.0.0"})
		conf = {"/static": 
					{"tools.staticdir.on": True,
					"tools.staticdir.dir": static_dir},
				"/results": 
					{"tools.staticdir.on": True,
					"tools.staticdir.dir": os.path.abspath(results_dir)}
				}
		# gather all of the data directories and make them available as static directories
		for k,v in data_dir.items():
			conf.update({"/" + k: 
							{"tools.staticdir.on": True,
							"tools.staticdir.dir": v}}) 

		app = Celery('grsnp')
		app.config_from_object('grsnp.celeryconfiguration')
		print "Checking for existing celery workers..."
		if app.control.inspect().ping() == None:
			if int(args['num_workers']) > 0:
				print "Starting Celery worker[s]..."
				fh = open("worker.log","w")
				script = ["celery","worker","-Q","long_runs,short_runs", '-c', str(args['num_workers']), "--app", "grsnp.worker_gr", "--loglevel", "INFO", "-n", "grsnp_LOCA1.%h",'-r',sett['run_files_dir'],'-d',args["data_dir"]]
				out = subprocess.Popen(script,stdout=fh,stderr=fh)
		else:
			workers = app.control.inspect().ping()
			pids = [str(app.control.inspect().stats()[j]['pid']) for j in workers.keys()]
			print "Celery workers already running. Pids:" + ",".join(pids) 
		cherrypy.config.update({'tools.sessions.timeout': 60})
		cherrypy.quickstart(WebUI(), "/", config=conf)

	else:
		print "WARNING: No port given. Server not started. Use --port flag to set port."


if __name__ == "__main__":	
	main()

