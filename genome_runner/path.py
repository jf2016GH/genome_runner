import os, fnmatch

from collections import defaultdict

def basename(k):
	return os.path.basename(os.path.splitext(k)[0])

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
