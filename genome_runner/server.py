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
import bedfilecreator as uscsreader
import query as grquery
from path import PathNode
from operator import itemgetter
from path import basename
import logging
from logging import FileHandler,StreamHandler
import json
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
		if os.listdir("results"):
			last_id = max(map(int, [basename(file) for file in os.listdir("results")]))
			while self.next_id() < last_id:
				pass

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
	def meta(self, tbl,organism="hg19"):
		try:
			trackdb = uscsreader.load_tabledata_dumpfiles("data/{}/trackDb".format(organism))
			html = trackdb[map(itemgetter('tableName'),trackdb).index(tbl)]['html']
		except Exception, e:
			print e
			return "<h3>ERROR: table not found</h3>"
		if html=='':
			return "<h3>(no data found).</h3>"
		else:
			return html

	@cherrypy.expose
	def query(self, bed_file=None,bed_data=None, background_file=None,background_data=None, niter=10, name="", score="", strand="", **kwargs):
		id = self.next_id()
		runset = {}
		cherrypy.response.timeout = 3600
		runset['filters'] = {"name": name,"score": score,"strand": strand}
		try:
			niter = int(niter)
		except ValueError:
			niter = 10
		runset["niter"] = niter

		try:
			jobname = kwargs["jobname"]
		except Exception, e:
			jobname = ""
			logger.error("id={}".format(id) + str(e))
		runset['jobname'] = jobname

			
		# load the FOI data
		bed_filename = ""
		data = ""
		try:
			f = os.path.join("uploads", str(id)+".bed")
			if not os.path.exists(f):
				with open(f, "wb") as out:
					if bed_file != None and bed_file.filename != "":
						bed_filename = bed_file.filename
						logger.info('Received uploaded FOI file (id={})'.format(id))
						while True:
							data = bed_file.file.read(8192)
							# TODO find empty lines
							#data = os.linesep.join([s for s in data.splitlines() if s ])

							# strips out new lines not compatible with bed tools
							#data = data.replace("\r\n","\n").replace("\r","\n")
							data = data.replace("\r","")
							if not data:
								break
							out.write(data)			
					elif bed_data!="":
						bed_filename = "Custom bed"
						logger.info('Received raw text  FOI data (id={})'.format(id))
						data = bed_data
						data = os.linesep.join([s for s in data.splitlines() if s])
						out.write(data)			
					else:
						return "upload a file please"
		except Exception, e:
			logger.error("id={}".format(id) + str(e))
			return "ERROR: upload a file please"
		runset["fois"] = bed_filename

		# load the background data if uploaded
		background_name = ""
		try:
			b = os.path.join("uploads",str(id) + ".background")
			if background_file != None and background_file.filename != "":
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
				background_name = "Custom bed"
				with open(b, "wb") as out:
					logger.info('Received raw text background data (id={})'.format(id))
					data = background_file
					data = os.linesep.join([s for s in data.splitlines() if s])
					out.write(data)
			else:
				background_name = "Genome"
				b = ""
		except Exception, e:
			logger.error("id={}".format(id) + str(e))
			return "ERROR: unable to upload background"
		runset['background'] = background_name

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


		for k,v in kwargs.items():
			# organism to use
			if "organism:" in v:
				organism = v.split(":")[-1]
			# which tests to run
			if "run:" in k and v=="on":
				print "change {}".format(k.split(":")[-1])
				run.append(k.split(":")[-1])

		runset['organism']= organism	
		runset['run'] = run

		# write the enrichment settings.
		path = os.path.join("results",str(id) + ".settings")
		settings = open(path,'wb')
		settings.write(json.dumps(runset))
		settings.close()


		# Reserve a spot in the results folder
		path = os.path.join("results", str(id))
		results = open(path,'wb') 
		results.close()

		# This starts the enrichment analysis in another OS process.
		# We know it is done when a file appears in the "results" directory
		# with the appropriate ID.
		p = Process(target=grquery.run_enrichments,
				args=(id,f,gfeatures,b,niter,name,score,strand,organism,run))
				
		p.start()
		raise cherrypy.HTTPRedirect("result?id=%d" % id)

	@cherrypy.expose
	def result(self, id):
		path = os.path.join("results", id)
		params = {}
		params["run_id"] = id
		if not os.path.exists(path):  #If file is empty...

			tmpl = lookup.get_template("enrichment_not_ready.html")
			return tmpl.render(id=id)
		if os.stat(path).st_size:
			with open(path) as strm:
				enrichments = cPickle.load(strm)
			enrichments.sort(key=attrgetter("p_value"))
		else:
			enrichments = []
		params["enrichments"] = enrichments
		tmpl = lookup.get_template("enrichment.html")
		progress = grquery.get_progress(id)
		# Loads the settings of the run if they exist
		if os.path.exists(path + ".settings"):
			settings = open(path + ".settings")
			s = json.loads(settings.read())
			settings.close()
			f = s['filters']
			# fills in default values for blank filter settings
			for key, value in s.iteritems():
				if value == "": s[key] = "None"
			for key, value in f.iteritems():
				if value == "": f[key] = "None"
			# merge all settings into one dictionary to pass to template
			params.update(s)
			params.update(f)
			print params
		else:
			params["background"] = "NA"
			params["name"] = "NA"
			params["niter"] = "NA"
			params["score"]= "NA"
			params["strand"] = "NA"
			params["foi"] = "NA"
			params["jobname"]= "NA"
		
		# Loads the progress file if it exists
		p = {"status":"","curprog":0,"progmax":0}
		if os.path.exists(path+".prog"):
			with open(path+".prog") as f:
				p = json.loads(f.read() )

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
	conf = {"/static": 
			{"tools.staticdir.on": True,
			"tools.staticdir.dir": static_dir}}
		
	cherrypy.quickstart(WebUI(), "/gr", config=conf)
