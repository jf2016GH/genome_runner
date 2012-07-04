import cherrypy, os, cgi, tempfile
import query as gquery
from mako.template import Template
from mako.lookup import TemplateLookup
lookup = TemplateLookup(directories=["templates"])

import fnmatch
from collections import defaultdict

from operator import attrgetter
from util import basename

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
					fs = ["file:"+os.path.join(base, f) for f 
						in list(fnmatch.filter(files, "*.bed"))]
					if fs:
						node[p] = fs

	def _li(self, k):
		if k.startswith("file:"):
			label = basename(k)
			name = k
		else:
			label = name = k
		return """\t<li><input name="%s" type="checkbox">
			<label>%s</label>\n""" % (name,label)

	def as_html(self, id=None):
		s = '<ul '
		if id:
			s += 'id="%s"' % str(id)
		s += '>'
		if not self.name == "Root":
			s += self._li(self.name)
		for k, child in sorted(self.items()):	
			if type(child) == list: #terminal
				s += "<ul>"
				for c in child:
					s += self._li(c)
				s += "</ul>"
			elif child:
				s += child.as_html()
		return s + "</ul>"

class WebUI(object):
	@cherrypy.expose
	def index(self):
		paths = PathNode()
		paths.name = "Root"
		paths.traverse("data")
		tmpl = lookup.get_template("index.html")
		return tmpl.render(paths=paths)

	@cherrypy.expose
	def query(self, bed_file=None, niter=10, **kwargs):
		try:
			niter = int(niter)
		except:
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

		gfeatures = [k[5:] for k,v in kwargs.items()
			if k.startswith("file:") and v=="on"]
		enrichments = []
		for gf in gfeatures:
			print f, gf
			e = gquery.enrichment(f,gf,n=niter)
			enrichments.append(e)
		enrichments.sort(key=attrgetter("p_value"))
		tmpl = lookup.get_template("enrichment.html")
		result = tmpl.render(enrichments=enrichments)
		os.unlink(f)
		return result

if __name__ == "__main__":
	cherrypy.server.max_request_body_size = 0
	static_dir = os.path.abspath(os.path.join(".", "static"))
	conf = {"/static": 
			{"tools.staticdir.on": True,
			"tools.staticdir.dir": static_dir}}
		
	cherrypy.quickstart(WebUI(), config=conf)
