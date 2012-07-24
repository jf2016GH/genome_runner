#!/usr/bin/env python
import sys, operator, os, numpy
from pybedtools import BedTool
from collections import namedtuple
import cPickle

from path import basename

# This class represents an Enrichment analysis result
# Lists of Enrichment objects are serialized to a Python Pickle
# file when an analysis is complete	
_Enrichment = namedtuple("Enrichment",
		["A","B","nA","nB","observed","expected","p_value","obsprox","expprox","pybed_p_value","pybed_expected"])
class Enrichment(_Enrichment):
	def category(self):
		if self.expected == 0 or self.p_value > 0.05:
			return "nonsig"
		elif self.observed < self.expected:
			return "under"
		else:
			return "over"

def make_filter(name, score, strand):
	def filter(interval):
		if name and name != interval.name:
			return False
		if score and score > int(interval.score):
			return False
		if strand and strand != interval.strand:
			return False
		return True
	return filter
# TODO make organism specific!
def enrichment(a, b,organism, name=None, score=None, strand=None, n=10):
	"""Perform enrichment analysis between two BED files.

	a - path to Feature of Interest BED file
	b - path to Genomic Feature BED file
	n - number of Monte-Carlo iterations
	"""
	organism = str(organism)
	print "Enrichment analysis organism is: {}".format(organism)
	flt = make_filter(name,score,strand)
	A = BedTool(str(a))
	B = BedTool(str(b)).filter(flt).saveas()
	nA = len(A)
	nB = len(B)
	if not nA or not nB:
		return Enrichment(a,basename(b),nA,nB,0,0,1,0,0,1,0)
	A.set_chromsizes(organism)
	B.set_chromsizes(organism)
	obs = len(A.intersect(B, u=True))
	# This is the Monte-Carlo step
	dist = [len(A.shuffle(genome=organism,chrom=True).intersect(B, u=True)) for i in range(n)]
	exp = numpy.mean(dist)
	p_value = len([x for x in dist if x > obs]) / float(len(dist))
	p_value = min(p_value, 1 - p_value)
	# expected calulated using pybed method
	pybeddist = A.randomintersection(B,iterations=n,shuffle_kwargs={'chrom': True})
	pybeddist = list(pybeddist)
	pybed_exp = numpy.mean(pybeddist)
	pybedp_value = len([x for x in pybeddist if x > obs]) / float(len(dist))
	pybedp_value = min(pybedp_value,1-pybedp_value)
	# expected calculated using pybed IntersectionMatrix
	#bed = [a]
	#im = IntersectionMatrix(beds,genome,iterations=n)
	#matrix = im.create_matrix(verbose=True)
	
	# proximety analysis
	# for expected

		#stores the means of the distances for the MC
	expall =[]
	for i in range(n):
		# run proximety analysis
		tmp = A.shuffle(genome=organism).closest(B,d=True)
		# get the distances
		for t in tmp:
			expall.append(t[-1])
	# calculate the overal expected distance
	expall.append(numpy.mean(numpy.array(expall,float)))
	# calculate the expected mean for all of the runs
	expprox = numpy.mean(numpy.array(expall,float))	
	# proximety analysis for observed
	tmp = A.closest(B,d=True)
	obsall = []
	for t in tmp:
		obsall.append(t[-1])
	obsprox = numpy.mean(numpy.array(obsall,float))

	return Enrichment(a, basename(b), nA, nB, obs, exp, p_value,obsprox,expprox,pybedp_value,pybed_exp)

def run_enrichments(id, f, gfeatures, niter, name, score, strand,organism):
	"""
	Run one FOI file (f) against multiple GFs, then 
	save the result to the "results" directory.
	"""
	enrichments = []
	for gf in gfeatures:
		e = enrichment(f,gf,organism,name,score,strand,niter)
		enrichments.append(e)
	path = os.path.join("results", str(id))
	with open(path, "w") as strm:
		cPickle.dump(enrichments, strm)

