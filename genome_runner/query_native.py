from ctypes import *

def tabulate(*args):
	return "\t".join(map(str, args))

class Enrichment(Structure):
	_fields_ = [
			("A", c_char_p),
			("B", c_char_p),
			("nA", c_int),
			("nB", c_int),
			("observed", c_int),
			("expected", c_double),
			("p_value", c_double)]

	def category(self):
		if self.p_value > 0.05 or self.expected == 0:
			return "notsig"
		ratio = self.observed / self.expected
		if ratio > 1:
			return "over"
		else:
			return "under"

	def __str__(self):
		return tabulate(self.nA, self.nB, self.p_value)

_libgr = cdll.LoadLibrary("./lib/libgenomerunner.so")
_libgr.enrichment.restype = Enrichment

def enrichment(A, B, n=10):
	return _libgr.enrichment("data/hg19.genome", A, B, n)

if __name__ == "__main__":
	g = "data/hg19.genome"
	a = "data/wgrna.bed"
	b = "data/rand.bed"
	print str(_libgr.enrichment(g,a,b,10))
