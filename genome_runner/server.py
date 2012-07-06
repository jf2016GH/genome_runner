import cherrypy, os, cgi, tempfile
import query as gquery
from mako.template import Template
from mako.lookup import TemplateLookup
lookup = TemplateLookup(directories=["templates"])

import fnmatch
from collections import defaultdict

from contextlib import closing
import sqlite3

from operator import attrgetter
from util import basename

#This class finds files recursively in the data/ directory
# so they can be represented as a treeview
class PathNode(defaultdict):
	def __init__(self):
		defaultdict.__init__(self, PathNode)
		self.files = []
	
	def traverse(self, base):
		for base, dirs, files in os.walk(base):
			prefix = base.split("/")[1:]
			node = self
			while prefix:
				p = prefix.pop(0)
				node = node[p]
				node.name = p
			node.files = ["file:"+os.path.join(base, f) for f 
				in list(fnmatch.filter(files, "*.bed"))]

	def _li(self, k):
		if k.startswith("file:"):
			label = basename(k)
			name = k
		else:
			label = name = k
		s = """\t<li><input name="%s" type="checkbox">
			<label>%s</label>\n""" % (name,label)
		if k.startswith("file:"):
			s = """<a target="_blank" href="meta?tbl=%s">%s</a>""" % \
				(label, s)
		return s

	def as_html(self, id=None):
		s = '<ul '
		if id:
			s += 'id="%s"' % str(id)
		s += '>'
		if not self.name == "Root":
			s += self._li(self.name)
		for f in sorted(self.files):	
			s += "<ul>" + self._li(f) + "</ul>"
		for _,child in sorted(self.items()):
			s += child.as_html()
		return s + "</ul>"

# Each function in this class is a web page
class WebUI(object):
	@cherrypy.expose
	def index(self):
		#Ultimately this HTML page will only need 
		# to be made once on startup
		#For testing this is nice because it 
		# instantly reflects new BED files added
		paths = PathNode()
		paths.name = "Root"
		paths.traverse("data")
		tmpl = lookup.get_template("index.html")
		return tmpl.render(paths=paths)

	@cherrypy.expose
	def meta(self, tbl):
		dbpath = "data/trackdb.db"
		with closing(sqlite3.connect(dbpath)) as conn:
			with closing(conn.cursor()) as c:
				c.execute("""select html from trackDb where
					tableName='%s'""" % tbl)
				try:
					s = str([h for (h,) in c][0])
					if not s:
						s = "<h3>(No data found).</h3>"
					return s
				except Exception, e:
					print e
					return "<h3>ERROR: table not found</h3>"

	@cherrypy.expose
	def query(self, bed_file=None, niter=10, **kwargs):
		# important! queries taking longer than 1h (3600s) will timeout
		# For such large queries it will be better to add a mechanism
		# To save query results to disk to be looked up later
		cherrypy.response.timeout = 3600
		try:
			niter = int(niter)
		except ValueError:
			niter = 10

		try:
			f = os.path.join("/tmp/", bed_file.filename)
			if not os.path.exists(f):

				with open(f, "w") as out:
					while True:
						data = bed_file.file.read(8192)
						if not data:
							break
						out.write(data)
		except:
			return "upload a file please"

		# "kwargs" (Keyword Arguments) stands for all the other
		# fields from the HTML form besides bed_file and niter
		# These other fields are all the tables whose boxes might
		# have been checked.
		gfeatures = [k[5:] for k,v in kwargs.items()
			if k.startswith("file:") and v=="on"]
		enrichments = []
		for gf in gfeatures:
			e = gquery.enrichment(f,gf,n=niter)
			enrichments.append(e)
		enrichments.sort(key=attrgetter("p_value"))
		tmpl = lookup.get_template("enrichment.html")
		result = tmpl.render(enrichments=enrichments)
		return result

if __name__ == "__main__":
	cherrypy.server.max_request_body_size = 0
	cherrypy.config.update({
		"server.socket_port":8081,
		"server.socket_host":"0.0.0.0"})
	static_dir = os.path.abspath(os.path.join(".", "static"))
	conf = {"/static": 
			{"tools.staticdir.on": True,
			"tools.staticdir.dir": static_dir}}
		
	cherrypy.quickstart(WebUI(), "/gr", config=conf)
