import cherrypy, os, cgi, tempfile
import query as gquery
from mako.template import Template
from mako.lookup import TemplateLookup
lookup = TemplateLookup(directories=["templates"])

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
		tmpl = lookup.get_template("index.html")
		return tmpl.render()

	@cherrypy.expose
	@cherrypy.tools.no_body_process()
	def query(self, bed_file=None, niter=10):
		try:
			niter = int(niter)
		except:
			niter = 10

		#Lots of voodoo here to handle file uploads.
		#See: http://tools.cherrypy.org/wiki/DirectToDiskFileUpload
		lcHDRS = {}
		for k,v in cherrypy.request.headers.iteritems():
			lcHDRS[k.lower()] = v
		fields = FieldStorageUL(fp=cherrypy.request.rfile,
				headers=lcHDRS,
				environ={"REQUEST_METHOD":"POST"},
				keep_blank_values=True)
		bed_file = fields["bed_file"]
		f = os.path.join("/tmp/", bed_file.filename)
		if not os.path.exists(f):
			os.link(bed_file.file.name, f)

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
