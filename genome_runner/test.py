import unpack_db
import os
import sys
from collections import Counter
from operator import itemgetter
#sys.stdout = open("myfile","w")

def printme():
	dbdir = os.path.join(os.environ["HOME"], "Work/hg19.db3")

	tracks = unpack_db.read_trackdb(dbdir)

	c = Counter(t.type for t in tracks if t.type[:3]=="bed")

	for k,v in sorted(c.items(), key=itemgetter(1)):
			print v, "\t", k
		
	#for k,v in sorted(c.items(), key=itemgetter(1)):
		#print v, "\t", k
		
"""
	for t in tracks:
		if "bed" in t.type[:3]:
			print t.type 

	for t in (x for x in tracks if x.type[:3] == "bed"):
		print(t.type)
		"""

def lose():
	print "you lose"

def win():
	print "you win"

def tie():
	print "it is a tie!"

outcomes = {"win":win,"lose":lose,"tie":tie}

def print_outcome(outcome):
	outcomes[outcome]()

def quarters(next_quarter=0.0):
	while True:
		yield next_quarter
		next_quarter += 0.25

result = []
for x in quarters():
	result.append(x)
	if x >= 1.0:
		break
	c = Counter(t for t in result if t>=0.2)
	print(c)
		
