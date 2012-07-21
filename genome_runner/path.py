import os, fnmatch

from collections import defaultdict
from collections import Set
def basename(k):
	return os.path.basename(os.path.splitext(k)[0])

#This class finds files recursively in the data/ directory
# so they can be represented as a treeview
class PathNode(defaultdict):

	# stores the organisms that are available for analysis.
	# each of the top level folders in the db directory are
	# considered to contain organism data.

	def __init__(self, organism="hg19"):
		defaultdict.__init__(self, PathNode)
		self.files = []
		self.organisms = []
	
	def traverse(self, base):
		for base, dirs, files in os.walk(base):
			prefix = base.split("/")[2:]
			node = self
			# creates a node for each subfolder i.e group,visibility
			while prefix:
				p = prefix.pop(0)
				if "." in p:
					continue
				node = node[p]
				node.name = p
			node.files = ["file:"+os.path.join(base, f) for f 
				in list(fnmatch.filter(files, ("*.gz")))]

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

	def as_html(self, id=None,):
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

	def org_as_html(self,id=None):
		s = "<select "
		if id:
			s += "id='%s'" % str(id)
		s += '>'
		for org in self.organisms:
			s += "<option>%s</option>" % str(org)
		return s + "</select>"

