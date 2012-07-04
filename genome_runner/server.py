import cherrypy, os, cgi, tempfile
import query as gquery
from mako.template import Template
from mako.lookup import TemplateLookup
lookup = TemplateLookup(directories=["templates"])

import fnmatch
from collections import defaultdict

class PathNode(defaultdict):
	def __init__(self):
		defaultdict.__init__(self, PathNode)
	
	def traverse(self, base):
		for base, dirs, files in os.walk(base):
			prefix = base.split("/")
			node = self
			while prefix:
				p = prefix.pop()
				if prefix:
					node = node[p]
					node.name = p
				else:
					fs = list(fnmatch.filter(files, "*.bed"))
					if fs:
						node[p] = fs

	def _li(self, k):
		return """\t<li><input name="%s" type="checkbox">
			<label>%s</label>\n""" % (k,k)

	def as_html(self, id=None):
		s = '<ul '
		if id:
			s += 'id="%s"' % str(id)
		s += '>'
		if not self.name == "Root":
			s += self._li(self.name)
		for k, child in self.items():			
			if type(child) == list: #terminal
				s += "<ul>"
				for c in child:
					s += self._li(c)
				s += "</ul>"
			elif child:
				s += child.as_html()
		return s + "</ul>"

class FieldStorageUL(cgi.FieldStorage):
	def make_file(self, binary=None):
		return tempfile.NamedTemporaryFile()

def no_body_process():
	cherrypy.request.process_request_body = False
cherrypy.tools.no_body_process = \
		cherrypy.Tool("before_request_body", no_body_process)


class WebUI(object):
	@cherrypy.expose
	def index(self):
		paths = PathNode()
		paths.name = "Root"
		paths.traverse("data")
		tmpl = lookup.get_template("index.html")
		return tmpl.render(paths=paths)

	@cherrypy.expose
	#@cherrypy.tools.no_body_process()
	def query(self, bed_file=None, niter=10, **kwargs):
		try:
			niter = int(niter)
		except:
			niter = 10

		f = os.path.join("/tmp/", bed_file.filename)
		with open(f, "w") as out:
			while True:
				data = bed_file.file.read(8192)
				if not data:
					break
				out.write(data)

		if os.path.exists(f):
			gfeatures = ["wgrna", "rand"]
			enrichments = []
			for gf in gfeatures:
				gf = "data/" + gf + ".bed"
				e = gquery.enrichment(f,gf,n=niter)
				enrichments.append(e)
			tmpl = lookup.get_template("enrichment.html")
			return tmpl.render(enrichments=enrichments)
		else:
			return "file not found"

if __name__ == "__main__":
	cherrypy.server.max_request_body_size = 0
	static_dir = os.path.abspath(os.path.join(".", "static"))
	conf = {"/static": 
			{"tools.staticdir.on": True,
			"tools.staticdir.dir": static_dir}}
		
	cherrypy.quickstart(WebUI(), config=conf)
