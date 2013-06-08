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

lookup = TemplateLookup(directories=["templates"])
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
	def __init__(self):
		self.next_id = itertools.count().next # This creates a threadsafe counter
		last_id = max(map(int, [file for file in os.walk('uploads').next()[1]]))
		while self.next_id() < last_id:
			pass
		print last_id

		self._index_html = {}

	@cherrypy.expose
	def index(self, organism="hg19"):
		if DEBUG_MODE or not organism in self._index_html:
			paths = PathNode(organism)
			paths.name = "Root"
			paths.organisms = self.get_org() 
			paths.traverse("data/{}".format(organism))
			tmpl = lookup.get_template("index.html")
			self._index_html[organism] = tmpl.render(paths=paths)
		return self._index_html[organism]

	def get_org(self):
		organisms = []
		files = os.listdir("./data")
		for f in files:
			if f.find(".") == -1:
				organisms.append(f)
		return organisms


	@cherrypy.expose
	def query(self, bed_file=None,bed_data=None, background_file=None,background_data=None, niter=10, name="", score="", strand="", **kwargs):
		id = self.next_id()
		print 'id: ', id
		upload_dir = os.path.join("uploads",str(id))
		os.mkdir(upload_dir)
		results_dir = os.path.join("results",str(id))
		os.mkdir(results_dir)
		fois = os.path.join(upload_dir,".fois") # contains a list of the paths to fois to run through the analysis
		gfs = os.path.join(upload_dir,".gfs") # contains a list of the paths to the gfs to run the fois against
		print id
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
						f = os.path.join(upload_dir, bed_filename)
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
							return "ERROR: An internal error has occured on the server."	
				# custom data entered	
				elif bed_data!="":
					f = os.path.join(upload_dir, "custom.bed")
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
		# load the background data if uploaded
		background_name = ""

		# "kwargs" (Keyword Arguments) stands for all the other
		# fields from the HTML form besides bed_file, niter, name, score, and strand
		# These other fields are all the tables whose boxes might
		# have been checked.
		# Thus with this way of doing things, it is not possible to have a genomicfeature
		# with one of these reserved names. 
		organism,run = "",[]
		gfeatures = [k[5:] for k,v in kwargs.items()
			if k.startswith("file:") and v=="on"]
		runset['gfs'] = gfeatures
		with open(gfs,"wb") as out_gfs:
			for g in gfeatures:
				out_gfs.write(g+"\n")

		for k,v in kwargs.items():
			# organism to use
			if "organism:" in v:
				organism = v.split(":")[-1]
			# which tests to run
			if "run:" in k and v=="on":
				print "change {}".format(k.split(":")[-1])
				run.append(k.split(":")[-1])

		runset['organism']= organism	
		runset['background'] = background_name
		
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
				#b_organ = {"hg19": "data/hg19/varRep/Tier0/snp137.gz",
				#		   "mm9": "data/mm9/varRep/Tier1/snp128.gz",
				#		   "mm8": "data/mm8/varRep/Tier1/snp126.gz"}
				background_name = "Genome"
				#b = b_organ[organism]
				b= ""
		except Exception, e:
			logger.error("id={}".format(id) + str(e))
			return "ERROR: unable to upload background"


		# write the enrichment settings.
		path = os.path.join(results_dir, ".settings")
		set_info = {"Jobname:": jobname,
					"Time:": strftime("%Y-%m-%d %H:%M:%S", gmtime()),
					"Background:": background_name,
					"Organism:": organism}
		with open(path, 'wb') as sett:
			for k,v in set_info.iteritems():
				sett.write(k+"\t"+v+"\n")

		# This starts the enrichment analysis in another OS process.
		# We know it is done when a file appears in the "results" directory
		# with the appropriate ID.
		p = Process(target=grquery.run_hypergeom,
				args=(fois,gfs,b,results_dir,runset['job_name'],True))				
		p.start()
		raise cherrypy.HTTPRedirect("result?id=%d" % id)

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
			print "progress_path: ", progress_path
			with open(progress_path) as f:
				p = json.loads(f.read() )

		# loads results from results file		
		detailed_path = os.path.join(path,"detailed.gr")
		if os.path.exists(detailed_path):
			with open(detailed_path) as f:
				params["detailed"] = f.read()

		# clustered matrix results loaded if they exist
		matrix_path = os.path.join(path,"matrix.gr")
		matrix_clust_path  = os.path.join(path,"matrix_clustered.gr")

		if os.path.exists(matrix_path):
			with open(matrix_path) as f:
				params["matrix"] = f.read().replace("\"","")
		if os.path.exists(matrix_clust_path):
			with open(matrix_clust_path) as f:
				params["matrix_data"] = f.read().replace("\"","")
				# d3 requires "gene_name" to be inserted into the first column
				tmp =  params["matrix_data"].split("\n")
				params["matrix_data"] = "\n".join(["\t".join(["gene_name",tmp[0]])]+tmp[1:])  
				params["matrix_data"] = params["matrix_data"].replace("\n","\\n")
		else: 
			params["matrix_data"] = "Heatmap will be available after the analysis is complete."

		params["log"] = "###Run Settings###\n"
		sett_path = os.path.join(path,".settings")
		if os.path.exists(sett_path):
			with open(sett_path) as f:
				params["log"] = params["log"]+ f.read()				
		params["log"] = params["log"] + "\n###Run Log###\n"
		debug_path = os.path.join(path,".log")
		if os.path.exists(debug_path):
			with open(debug_path) as f:
				params["log"] = params["log"] + f.read()

		# check if run files ready for download
		zip_path = os.path.join(path,"GR_Runfiles_{}.zip".format([y.split("\t")[1] for y in open(sett_path).read().split("\n") if "Jobname:" in y][0]))
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
	conf = {"/static": 
			{"tools.staticdir.on": True,
			"tools.staticdir.dir": static_dir},
			"/results": 
			{"tools.staticdir.on": True,
			"tools.staticdir.dir": results}
			}
		
	cherrypy.quickstart(WebUI(), "/gr", config=conf)
