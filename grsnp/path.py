import os, fnmatch,json

from collections import defaultdict
from collections import Set
from collections import namedtuple
import re	



#This class finds files recursively in the data/ directory
# so they can be represented as a treeview
class PathNode(defaultdict):

	# stores the organisms that are available for analysis.
	# each of the top level folders in the db directory are
	# considered to contain organism data.

	static_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),"frontend","static")

	def __init__(self):
		defaultdict.__init__(self, PathNode)
		self.files = []
		self.organisms = []

	# for the autocomplete text box
	def traverse(self, base):
		''' Reads the data directory for GenomicFeatures
		'''
		blacklist = []
		blacklist_path = os.path.join(base,"blacklist.txt")
		gfs_path = os.path.join(base,"gfs.php")

		if os.path.exists(blacklist_path):
			with open(blacklist_path) as f:
				blacklist = [line.strip() for i,line in enumerate(f)]
		# used to generate a json list of gfs
		gfs = []
		int_data = len(base.split("/"))
		for base, dirs, files in os.walk(base):

			prefix = base.split("/")[int_data:]
			node = self
			# creates a node for each subfolder i.e group,visibility
			while prefix:
				p = prefix.pop(0)
				if "." in p:
					continue
				node = node[p]
				node.name = p
			node.files = ["file:"+os.path.join(base, f) for f 
				in files if f.endswith(('.gz', '.bb')) and base_name(f) not in blacklist]

			# used for the auto-complete text box
			for f in node.files:
				if base_name(f) not in blacklist:
					 gfs.append({"caption": str(base_name(os.path.splitext(f)[0])),"value":f})
		f = open(gfs_path,"wb")
		f.write(json.dumps(gfs))
		f.close()

    
	def _li(self, k):
		if k.startswith("file:"):
			label = base_name(k)
			name = k
		else:
			label = name = k
		s = """\t<li><input name="%s" type="checkbox">
			<label>%s</label>\n""" % (name,label)
		if k.startswith("file:"):
			org= k.split(os.sep)[k.split(os.sep).index("grsnp_db")+1]
			s = """<a target="_blank" href="meta?tbl=%s&organism=%s">%s</a>""" % \
				(label,org, s)

		s = """\t<li><input name="%s" type="checkbox">
			<label>%s</label>\n""" % (name,label)
		if k.startswith("file:"):
			org= k.split(os.sep)[k.split(os.sep).index("grsnp_db")+1]
			s = """<a target="_blank" href="meta?tbl=%s&organism=%s">%s</a>""" % \
				(label,org, s)
		return s

	def _treeview_html(self, id=None,):
		s = '<ul '
		if id:
			s += 'id="%s"' % str(id)
		s += '>'
		if not self.name == "Root":
			s += self._li(self.name)
		for f in sorted(self.files):	
			# ensure that trackDb is not listed
			if "trackDb" not in f and "chromInfo" not in f:
				s += "<ul>" + self._li(f) + "</ul>"
		for _,child in sorted(self.items()):
			s += child._treeview_html()		
		return s + "</ul>"


	def write_treeview_html(self,base,organism):
		html_path = os.path.join(base,organism,"treeview.html")
		html = self._treeview_html(id="ucsc") # generate the treeview html code
		with open(html_path,"wb") as writer:
			writer.write(html)

	# for the organism combobox
	def org_as_html(self,id=None):

		s = "<select name='organism' style=\"margin-left: 10px\" "
		if id:
			s += "id='%s'" % str (id)
		s += '>'
		for org in self.organisms:
			s += "<option value='organism:{}'>{}</option>".format(org,org)
		return s + "</select>"


	def get_custom_fois(self,organism,custom_dir):
		'''Get all of the sample FOI SNPs files.
		'''
		demo_dir = os.path.join(custom_dir,"fois",organism)
		html = """<button type="button" id="demo_fois_none" onclick="enable_foi_uploads()" style="margin-top: 12px"  class="btn btn-primary active" title="" >None</button>\n"""
		if not os.path.exists(demo_dir):
			return html
		for snp_dir in [ os.path.join(demo_dir,f) for f in os.listdir(demo_dir) if os.path.isdir(os.path.join(demo_dir,f))]:
			tooltip = "Includes the following files:\n"
			for s in [os.path.join(snp_dir,f) for f in os.listdir(snp_dir) if os.path.isfile(os.path.join(snp_dir,f)) and not f.endswith(".tbi")]:
				tooltip += "\t"+base_name(s) + "\n" 
			html = html + """<button type="button" onclick="clear_foi_uploads()" style="margin-top: 12px"  class="btn btn-primary" data-toggle-value="{}" title="{}" >{}</button>\n""".format(snp_dir,tooltip,base_name(snp_dir))
		return html

	def get_backgrounds_combo(self,organism,custom_dir):
		''' Generates the html code for the combo box containing the 
			default organism backgrounds.
		'''

		html = """<select name="default_background" style="margin-left: 5px; margin-top: 9px" id="default_background">"""
		background_dir = os.path.join(custom_dir,"backgrounds",organism)
		if not os.path.exists(background_dir):
			return html + "</select>"
		for bk in [ f for f in sorted(os.listdir(background_dir)) if os.path.isfile(os.path.join(background_dir,f)) and not f.endswith(".tbi")]:
			tmp = os.path.join(background_dir,bk)			
			html = html + "<option value='{}'>{}</option>".format(tmp,re.sub("^\d+",'',tmp.split("/")[-1].split(".")[0]))
		html  = html + "</select>"
		return html

	def get_custom_gfs(self,organism,custom_dir):
		demo_dir = os.path.join(custom_dir,"gfs",organism)
		html = ""
		if not os.path.exists(demo_dir):
			return ""
		for gfs_dir in [ os.path.join(demo_dir,f) for f in os.listdir(demo_dir) if os.path.isdir(os.path.join(demo_dir,f))]:
			tooltip = "Includes the following files:\n"
			for s in [os.path.join(gfs_dir,f) for f in os.listdir(gfs_dir) if os.path.isfile(os.path.join(gfs_dir,f)) and not f.endswith(".tbi")]:
				tooltip += "\t"+base_name(s) + "\n" 
			html += """<input type="checkbox" style="font-size:120%;"  name="grouprun:{}" style="margin: 10px">{}</input>
						<img class="helptooltip" title="{}" style="position: relative;top: 6px;" width="25" height="25" src="static/images/help-icon.png" alt="help">""".format(gfs_dir,base_name(gfs_dir),tooltip)
		return html

def base_name(k):
    return os.path.basename(k).split(".")[0]