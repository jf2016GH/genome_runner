#!/usr/bin/env python
import sys, operator, os, numpy
from pybedtools import BedTool
import pybedtools
from collections import namedtuple
import cPickle



from path import basename

# This class represents an Enrichment analysis result
# Lists of Enrichment objects are serialized to a Python Pickle
# file when an analysis is complete	
_Enrichment = namedtuple("Enrichment",
		["A","B","nA","nB","observed","expected","p_value","obsprox","expprox","pybed_p_value","pybed_expected","jaccard_observed","jaccard_p_value","jaccard_expected"])
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
	print A.count()
	print B.count()
	nA = len(A)
	nB = len(B)
	if not nA or not nB:
		return Enrichment(a,basename(b),nA,nB,0,0,1,0,0,1,0,0,1,0)
	A.set_chromsizes(organism)
	B.set_chromsizes(organism)
	obs = len(A.intersect(B, u=True))
	# This is the Monte-Carlo step
	print "RUNNING MONTE CARLO"
	dist = [len(A.shuffle(genome=organism,chrom=True).intersect(B, u=True)) for i in range(n)]
	exp = numpy.mean(dist)
	# gave p_value a value here so that it doesn't go out of scope, is this needed?
	p_value = 'NA'
	if exp == obs or (exp == 0 and obs == 0):
		p_value =1
	else:
		p_value = len([x for x in dist if x > obs]) / float(len(dist))
		p_value = min(p_value, 1 - p_value)
	
	print "MC pvalue: {}".format(p_value)
	# expected caluclated using pybed method
	print "RUNNING RANDOM INTERSECTIONS"
	pybeddist = A.randomintersection(B,iterations=n,shuffle_kwargs={'chrom': True})
	pybeddist = list(pybeddist)
	print pybeddist
	print "calculating mean"
	pybed_exp = numpy.mean(pybeddist)
	pybedp_value = 'NA'
	if pybed_exp == obs or (pybed_exp == 0 and obs == 0):
		pybedp_value =1
	else:
		print "calculating pvalue"
		pybedp_value = len([x for x in pybeddist if x > obs]) / float(len(pybeddist))
		pybedp_value = min(pybedp_value,1-pybedp_value)
	# epected calculated using jaccard method
	chrom = pybedtools.get_chromsizes_from_ucsc(organism)
	chrom_fn = pybedtools.chromsizes_to_file(chrom)
	A2 = A.cut([0,1,2])
	B2 = B.cut([0,1,2])
	print "RUNNING JACCARD"
	resjaccard = A2.naive_jaccard(B2,genome_fn=chrom_fn,iterations=n,shuffle_kwargs={'chrom':True})
	jaccard_dist = resjaccard[1]
	jaccard_obs = resjaccard[0]
	jaccard_exp = numpy.mean(resjaccard[1])
	jaccardp_value = 'NA'
	if jaccard_exp == jaccard_obs or (jaccard_exp  == 0 and jaccard_obs == 0):
		jaccardp_value =1
	else:
		jaccardp_value = len([x for x in jaccard_dist if x > obs]) / float(len(jaccard_dist))
		jaccardp_value = min(pybedp_value,1-pybedp_value)
	print "JACCARD PVALUE: {}".format(jaccardp_value)

	# run proximety analysis
	if True:
		print "RUNNING PROXIMITY"
			#stores the means of the distances for the MC
		expall =[]
		for i in range(n):
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
	else:
		print "SKIPPING PROXIMITY"
		obsprox = -1
		expprox = -1
	print "FINISHED"

	return Enrichment(a, basename(b), nA, nB, obs, exp, p_value,obsprox,expprox,pybedp_value,pybed_exp,jaccard_obs,jaccardp_value,jaccard_exp)

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

