import cherrypy, os, cgi, tempfile, sys, itertools
from mako.template import Template
from mako.lookup import TemplateLookup

from contextlib import closing
import sqlite3
from operator import attrgetter
from multiprocessing import Process
import cPickle
import bedfilecreator as uscsreader
import query
from path import PathNode
from operator import itemgetter

lookup = TemplateLookup(directories=["templates"])
DEBUG_MODE = True

# Each function in this class is a web page
class WebUI(object):
	def __init__(self):
		self.next_id = itertools.count().next # This creates a threadsafe counter
		if os.listdir("results"):
			last_id = max(map(int, os.listdir("results")))
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
	def query(self, bed_file=None,bed_data=None, niter=10, name="", score="", strand="", **kwargs):
		cherrypy.response.timeout = 3600
		try:
			niter = int(niter)
		except ValueError:
			niter = 10

			
		data = ""
		try:
			id = self.next_id()
			f = os.path.join("uploads", str(id)+".bed")
			if not os.path.exists(f):
				with open(f, "wb") as out:
					if bed_file != None and bed_file.filename != "":
						while True:
							data = bed_file.file.read(8192)
							data = data.strip("\r")
							data = os.linesep.join([s for s in data.splitlines() if s])
							if not data:
								break
							out.write(data)			
					elif bed_data!="":
						data = bed_data
						data = os.linesep.join([s for s in data.splitlines() if s])
						out.write(data)			
					else:
						return "upload a file please"

		except Exception, e:
			print e
			return "ERROR: upload a file please"
		print "### LOADED DATA####"
		# "kwargs" (Keyword Arguments) stands for all the other
		# fields from the HTML form besides bed_file, niter, name, score, and strand
		# These other fields are all the tables whose boxes might
		# have been checked.
		# Thus with this way of doing things, it is not possible to have a genomicfeature
		# with one of these reserved names. 
		organism = ""
		gfeatures = [k[5:] for k,v in kwargs.items()
			if k.startswith("file:") and v=="on"]
		for k,v in kwargs.items():
			if v.startswith("organism:"):
				print "organism is: " + k 
				organism = v.split(":")[-1]
				print organism

		# Reserve a spot in the results folder
		path = os.path.join("results", str(id))
		results = open(path,'wb') 
		results.close()
		
		# This starts the enrichment analysis in another OS process.
		# We know it is done when a file appears in the "results" directory
		# with the appropriate ID.
		p = Process(target=query.run_enrichments,
				args=(id,f,gfeatures,niter,name,score,strand,organism))
		p.start()
		raise cherrypy.HTTPRedirect("result?id=%d" % id)

	@cherrypy.expose
	def result(self, id):
		path = os.path.join("results", id)
		if not os.path.exists(path):
			tmpl = lookup.get_template("enrichment_not_ready.html")
			return tmpl.render(id=id)
		with open(path) as strm:
			enrichments = cPickle.load(strm)
		enrichments.sort(key=attrgetter("p_value"))
		tmpl = lookup.get_template("enrichment.html")
		return tmpl.render(enrichments=enrichments)


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
