import cherrypy, os, cgi, tempfile, sys, itertools
from mako.template import Template
from mako.lookup import TemplateLookup
from mako.exceptions import RichTraceback as MakoTraceback
from contextlib import closing
import sqlite3
import re
from operator import attrgetter
from multiprocessing import Process
import cPickle
from path import PathNode
from operator import itemgetter
from path import basename
import logging
from logging import FileHandler,StreamHandler
import json
import pdb
import hypergeom3 as grquery
from time import gmtime, strftime
import bedfilecreator
import simplejson
import string
import random
import traceback
import bedfilecreator as uscsreader

lookup = TemplateLookup(directories=["templates"])


sett = {"dir_data": "data",
		"default_organism": "hg19"
		}

# add default background paths here
default_backgrounds_paths = {"hg19": [os.path.join(sett["dir_data"],"hg19/genes/Tier1/all_diseases2.gz"),
									  os.path.join(sett["dir_data"],"hg19/genes/Tier2/all_diseases1.gz")],
							  "mm9": "data/mm9/varRep/Tier1/snp128.gz",
							  "mm8": "data/mm8/varRep/Tier1/snp126.gz"}
DEBUG_MODE = True

logger = logging.getLogger('genomerunner.server')
hdlr = logging.FileHandler('genomerunner_server.log')
hdlr_std = StreamHandler()
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.addHandler(hdlr_std)
logger.setLevel(logging.INFO)


# Each function in this class is a web page 
class WebUI(object):
	def __init__(self, default_organism="hg19"):
		self._index_html = {}

	@cherrypy.expose
	def index(self, organism=sett["default_organism"]):
		if DEBUG_MODE or not organism in self._index_html:
			paths = PathNode()
			paths.name = "Root"
			paths.organisms = self.get_org() 
			paths.traverse(os.path.join(sett["dir_data"],organism))
			tmpl = lookup.get_template("index.html")
			# Load default backgrounds


			self._index_html[organism] = tmpl.render(paths=paths,default_background=self.get_backgrounds_combo(organism))
		return self._index_html[organism]


	def get_org(self):
		organisms = []
		files = os.listdir(sett["dir_data"])
		for f in files:
			if f.find(".") == -1:
				organisms.append(f)
		return organisms

	def get_backgrounds_combo(self,organism):
		''' Generates the html code for the combo box containing the 
			default organism backgrounds.
		'''

		html = """<select name="default_background" style="margin-left: 5px; margin-top: 9px" id="default_background">"""
		for bk in default_backgrounds_paths[organism]:			
			html = html + "<option value='{}'>{}</option>".format(bk,bk.split("/")[-1].split(".")[0])
		html  = html + "</select>"
		return html

	@cherrypy.expose
	def query(self, bed_file=None,bed_data=None, background_file=None,background_data=None, 
				genomicfeature_file=None, niter=10, name="", score="", strand="",run_annotation=None, default_background = "",**kwargs):
		# Assign a random id
		id = ''.join(random.choice(string.lowercase+string.digits) for _ in range(32))
		while (os.path.exists(os.path.join("uploads",id))):
			id = ''.join(random.choice(string.lowercase+string.digits) for _ in range(32))

		upload_dir = os.path.join("uploads",str(id))
		os.mkdir(upload_dir)
		os.mkdir(os.path.join(upload_dir,"fois"))
		os.mkdir(os.path.join(upload_dir,"gfs"))
		results_dir = os.path.join("results",str(id))
		os.mkdir(results_dir)
		fois = os.path.join(upload_dir,".fois") # contains a list of the paths to fois to run through the analysis
		gfs = os.path.join(upload_dir,".gfs") # contains a list of the paths to the gfs to run the fois against
		runset = {}
		cherrypy.response.timeout = 3600

		try:
			jobname = kwargs["jobname"]
		except Exception, e:
			jobname = ""
			logger.error("id={}".format(id) + str(e))
		runset['job_name'] = jobname
		runset['time'] = strftime("%Y-%m-%d %H:%M:%S", gmtime())

			
		# load the FOI data
		bed_filename = ""

		data = ""
		try:
			with open(fois,"wb") as out_fois:
				# bed files uploaded
				if bed_file:
					if not isinstance(bed_file,(list)): bed_file = [bed_file] # makes a list if only one file uploaded
					for b in bed_file:
						bed_filename = b.filename
						f = os.path.join(upload_dir, "fois",bed_filename)
						if not os.path.exists(f):
							with open(f, "wb") as out:
								if b != None and b.filename != "":
									logger.info("Received uploaded FOI file (name={}, id={})".format(bed_filename, id))
									while True:
										data = b.file.read(8192)
										# TODO find empty lines
										#data = os.linesep.join([s for s in data.splitlines() if s ])

										# strips out new lines not compatible with bed tools
										data = data.replace("\r","")
										if not data:
											break
										out.write(data)			
									out_fois.write(f+"\n")
						else:
							logger.error("id={} Upload file already exists at {}".format(id,f))
							print "id={} Upload file already exists at {}".format(id,f)
				# custom data entered	
				elif bed_data!="":
					f = os.path.join(upload_dir,"fois", "custom.bed")
					with open(f, "wb") as out:
						bed_filename = "custom.bed"
						logger.info('Received raw text  FOI data (id={})'.format(id))
						data = bed_data
						data = os.linesep.join([s for s in data.splitlines() if s])
						out.write(data)		
					out_fois.write(f+"\n")	
				else:
					return "upload a file please"

		except Exception, e:
			logger.error("id={}".format(id) + str(e))
			return "ERROR: upload a file please"
		runset["fois"] = bed_filename

		# uploads custom genomic features
		try:
			with open(gfs,"wb") as out_gfs:
				# bed files uploaded
				if genomicfeature_file:
					if not isinstance(genomicfeature_file,(list)): genomicfeature_file = [genomicfeature_file] # makes a list if only one file uploaded
					for b in genomicfeature_file:
						gfbed_filename = b.filename
						f = os.path.join(upload_dir, "gfs", gfbed_filename)
						if not os.path.exists(f):
							with open(f, "wb") as out:
								if b != None and b.filename != "":
									logger.info("Received uploaded GF file (name={}, id={})".format(gfbed_filename, id))
									while True:
										data = b.file.read(8192)
										# TODO find empty lines
										#data = os.linesep.join([s for s in data.splitlines() if s ])

										# strips out new lines not compatible with bed tools
										data = data.replace("\r","")
										if not data:
											break
										out.write(data)			
									out_gfs.write(f+"\n")
						else:
							logger.error("id={} Uploaded GF file already exists at {}".format(id,f))
		except Exception, e:
			logger.error("id={}".format(id) + str(e))
			return "ERROR: Unable to process custom Genome annotation feature"

		# "kwargs" (Keyword Arguments) stands for all the other
		# fields from the HTML form besides bed_file, niter, name, score, and strand
		# These other fields are all the tables whose boxes might
		# have been checked.
		# Thus with this way of doing things, it is not possible to have a genomicfeature
		# with one of these reserved names. 
		organism,run = "",[]
		gfeatures = [k[5:] for k,v in kwargs.items()
			if k.startswith("file:") and v=="on"]
		with open(gfs,"a") as out_gfs:
			for g in gfeatures:
				out_gfs.write(g+"\n")

		for k,v in kwargs.items():
			# organism to use
			if "organism:" in v:
				organism = v.split(":")[-1]
			# which tests to run
			if "run:" in k and v=="on":
				run.append(k.split(":")[-1])
			# append custom list to be run
			if "grouprun:" in k and v == "on":
				gp_gfs = open( k.split(":")[-1],"rb").read()
				with open(gfs,"a") as out_gfs:
					out_gfs.write(gp_gfs)

		runset['organism']= organism	
		
		# load the background data if uploaded
		background_name = ""
		try:
			if background_file != None and background_file.filename != "":
				b = os.path.join(upload_dir,background_file.filename+".background")
				logger.info('Received uploaded background file (id={})'.format(id))
				background_name = background_file.filename
				with open(b, "wb") as out:
					while True:
						data = background_file.file.read(8192)
						# TODO find empty lines
						#data = os.linesep.join([s for s in data.splitlines() if not s.isspace() ])

						# strips out new lines not compatible with bed tools
						data = data.replace("\r","")
						if not data:
							break
						out.write(data)			
			elif background_data != None and background_data != "":
				b = os.path.join(upload_dir,"custom.background")
				background_name = "custom.bed"
				with open(b, "wb") as out:
					logger.info('Received raw text background data (id={})'.format(id))
					data = background_file
					data = os.linesep.join([s for s in data.splitlines() if s])
					out.write(data)
			else:				
				background_name = default_background.split("/")[-1].split(".")[0] 
				b = default_background
		except Exception, e:
			logger.error("id={}".format(id) + str(e))
			return "ERROR: unable to upload background"
		runset['background'] = background_name


		# write the enrichment settings.
		path = os.path.join(results_dir, ".settings")
		set_info = {"Jobname:": jobname,
					"Time:": strftime("%Y-%m-%d %H:%M:%S", gmtime()),
					"Background:": background_name,
					"Organism:": organism,
					"Run Annotation:": str(bool(run_annotation))}

		with open(path, 'wb') as sett:
			for k,v in set_info.iteritems():
				sett.write(k+"\t"+v+"\n")

		# This starts the enrichment analysis in another OS process.
		# We know it is done when a file appears in the "results" directory
		# with the appropriate ID.
		p = Process(target=grquery.run_hypergeom,
				args=(fois,gfs,b,results_dir,runset['job_name'],True))				
		p.start()
		raise cherrypy.HTTPRedirect("result?id=%s" % id)

	@cherrypy.expose
	def result(self, id):
		path = os.path.join("results", id)
		params = {}
		params["run_id"] = id
		params["detailed"] = "Results not yet available"
		params["matrix"] = "Results not yet available"
		if not os.path.exists(path):  #If file is empty...
			tmpl = lookup.get_template("enrichment_not_ready.html")
			return tmpl.render(id=id)
		tmpl = lookup.get_template("results.html")

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

		
		foi_names_path = os.path.join(os.path.join("uploads", id),".fois")
		if os.path.exists(foi_names_path):
			with open(foi_names_path) as f:
				params["fois"] = [basename(x).split(".")[0] for x in f.read().split("\n") if x != ""]
		else:
			params["fois"] = ""
		# check if run files ready for download
		zip_path = os.path.join(path,"GR_Runfiles_{}.tar.gz".format([y.split("\t")[1] for y in open(sett_path).read().split("\n") if "Jobname:" in y][0]))
		if os.path.exists(zip_path):
			params["zipfile"] = zip_path
		else:
			params["zipfile"] = ""

		params.update(p)
		try:
			rend_template = tmpl.render(**params)
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
	def get_heatmaps(self, run_id, organism):
		"""	Returns clustered and PCC matrix if they exist.
		'organism': is used to load detailed labels for the GFs.
		"""
		cherrypy.response.headers['Content-Type'] = 'application/json'
		trackdb = bedfilecreator.load_tabledata_dumpfiles(os.path.join("data",organism,"trackDb"))
		results = {}
		path = os.path.join("results", run_id)
		matrix_path = os.path.join(path,"matrix.txt")
		# Load clustered matrix
		matrix_clust_path  = os.path.join(path,"clustered.txt")
		if os.path.exists(matrix_path):
			with open(matrix_path) as f:
				results["matrix"] = f.read().replace("\"","")
		if os.path.exists(matrix_clust_path):
			with open(matrix_clust_path) as f:
				d = f.read()
				results["matrix_data"] = d.replace("\"","")
				# d3 requires "gene_name" to be inserted into the first column
				tmp_data =  results["matrix_data"].split("\n")
				tmp = tmp_data[:] # this copy is passed onto the results page
				# insert the description column if it does not exist
				tmp_matrix = [x.split("\t") for x in tmp]
				if tmp_matrix[0][-1] != "Genomic Feature Description":
					tmp_matrix[0] += ["Genomic Feature Description"]
					for i in range(1,len(tmp)):
						description = [x["longLabel"] for x in trackdb if x["tableName"] == tmp_matrix[i][0]]
						if len(description) is not 0: description = description[0]
						else: description = ""
						tmp_matrix[i] += [description]						

			results["matrix_data"] = "\n".join(["\t".join(["gene_name",tmp[0]])]+tmp[1:])  
			results["matrix_data"] = results["matrix_data"]
			results["matrix_data_gf_description"] = "\t".join([x[-1] for x in tmp_matrix[1:]])
		else: 
			results["matrix_data"] = "Heatmap will be available after the analysis is complete."
			results["matrix_data_gf_description"] = ""

		# Pearson's matrix results
		matrix_cor_path = os.path.join(path,"pcc_matrix.txt")
		if os.path.exists(matrix_cor_path):
			with open(matrix_cor_path) as f:
				d = f.read()
				results["matrix_cor_data"] =  d.replace("\"","")
				results["matrix_cor"] = d.replace("\"","")
				# d3 requires "gene_name" to be inserted into the first column
				tmp =  results["matrix_cor"].split("\n")
				results["matrix_cor"] = "\n".join(["\t".join(["gene_name",tmp[0]])]+tmp[1:])  
				results["matrix_cor"] = results["matrix_cor"]				

		else:
			results["matrix_cor_data"] = ""
			results["matrix_cor"] = "Heatmap wil-l be available after the analysis is complete."
		pvalue_path = os.path.join(path,"pcc_matrix_pvalue.txt")
		if os.path.exists(pvalue_path):
			with open(pvalue_path) as f:
				d = f.read()
				results["matrix_cor_pvalues"] = d.replace("\"","")
				# d3 requires "gene_name" to be inserted into the first column
				tmp =  results["matrix_cor_pvalues"].split("\n")
				results["matrix_cor_pvalues"] = "\n".join(["\t".join(["gene_name",tmp[0]])]+tmp[1:])   
				results["matrix_cor_pvalues"] = results["matrix_cor_pvalues"]			
		else: 
			results["matrix_cor_pvalues"] = ""
		return simplejson.dumps(results)


	@cherrypy.expose
	def get_annotation(self,run_id,foi_name):
		annotation_path = os.path.join(os.path.join("results",run_id,foi_name + ".txt"))
		results = []
		if os.path.exists(annotation_path):
			with open(annotation_path) as f:
				# skip the comment lines
				cols = f.readline().rstrip()
				while cols[0] == "#":
					cols = f.readline().rstrip()
				cols = cols.split("\t")	
				results.append(cols)			
				for foi in f:
					if foi.strip() != "":
						results.append(foi.rstrip().split("\t"))
		return simplejson.dumps(results)

	@cherrypy.expose
	def meta(self, tbl,organism="hg19"):
		try:
			trackdb = uscsreader.load_tabledata_dumpfiles("data/{}/trackDb".format(organism))
			html = trackdb[map(itemgetter('tableName'),trackdb).index(tbl)]['html']
		except Exception, e:
			print traceback.format_exc()
			return "<h3>(No data found for {}.)</h3>".format(tbl)
		if html=='':
			return "<h3>(No data found for {}.)</h3>".format(tbl)
		else:
			return html

	@cherrypy.expose
	def get_detailed(self,run_id):
		""" loads results from detailed results file
		"""
		detailed_path,results = os.path.join("results", run_id,"detailed.txt"),{"detailed": ""}		 
		if os.path.exists(detailed_path):
			with open(detailed_path) as f:
				results["detailed"] = f.read()
		return simplejson.dumps(results)

	@cherrypy.expose
	def get_progress(self, run_id):
		# Loads the progress file if it exists
		p = {"status":"","curprog":0,"progmax":0}
		progress_path = os.path.join(os.path.join("results", run_id),".prog")
		if os.path.exists(progress_path):
			with open(progress_path) as f:
				p = json.loads(f.read() )
		return simplejson.dumps(p)

	@cherrypy.expose
	def get_log(sefl,run_id):
		results = {"log": ""}
		log_path = os.path.join(os.path.join("results", run_id),"log.txt")
		print log_path
		if os.path.exists(log_path):
			with open(log_path) as f:
				results["log"] = f.read()
				print "log ",results["log"]
		return simplejson.dumps(results)

	@cherrypy.expose
	def enrichment_log(self, id):
		with open(os.path.join("results",id+".log")) as sr:
			x = sr.read()
			return "<p>{}</p>".format(x.replace("\n","<br/>"))

	@cherrypy.expose
	def cite(self):
		tmpl = lookup.get_template("cite.html")
		return tmpl.render()

	@cherrypy.expose
	def news(self):
		tmpl = lookup.get_template("news.html")
		return tmpl.render()

	@cherrypy.expose
	def overview(self):
		tmpl = lookup.get_template("overview.html")
		return tmpl.render()

	@cherrypy.expose
	def demo(self):
		tmpl = lookup.get_template("demo.html")
		return tmpl.render()

if __name__ == "__main__":
	if not os.path.exists("results"):
		os.mkdir("results")
	if not os.path.exists("uploads"):
		os.mkdir("uploads")

	if len(sys.argv) >= 2:
		port = int(sys.argv[1])
	else:
		port = 8081
	cherrypy.server.max_request_body_size = 0
	cherrypy.config.update({
		"server.socket_port":port,
		"server.socket_host":"0.0.0.0"})
	static_dir = os.path.abspath(os.path.join(".", "static"))
	results = os.path.abspath(os.path.join(".","results"))
	media = os.path.abspath(os.path.join(".","media"))
	conf = {"/static": 
				{"tools.staticdir.on": True,
				"tools.staticdir.dir": static_dir},
			"/results": 
				{"tools.staticdir.on": True,
				"tools.staticdir.dir": results},
			"/media":
				{"tools.staticdir.on": True,
				"tools.staticdir.dir": media},
			}
		
	cherrypy.quickstart(WebUI(), "/gr", config=conf)
