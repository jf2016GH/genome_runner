#!/usr/bin/env python
import sys, operator, os, numpy
from pybedtools import BedTool
from collections import namedtuple
from util import basename

Enrichment = namedtuple("Enrichment",
		["A","B","nA","nB","observed","expected","p_value"])
def enrichment_category(self):
	if self.p_value > 0.05:
		return "nonsig"
	elif self.observed < self.expected:
		return "under"
	else:
		return "over"
Enrichment.category = enrichment_category

def monte_carlo(ref, query, n=10):
	result = []
	for i in range(n):
		shuffled = query.shuffle(g="data/hg19.genome")
		overlaps = shuffled.intersect(ref,u=True)
		result.append(len(overlaps))
	return result

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

def enrichment(a, b, name=None, score=None, strand=None, n=10):
	"""Perform enrichment analysis between two FeatureSets.

	a - Feature of Interest BedTool object
	b - Genomic Feature BedTool object
	"""
	flt = make_filter(name,score,strand)
	A = BedTool(str(a)).filter(flt).saveas()
	B = BedTool(str(b)).filter(flt).saveas()
	nA = len(A)
	nB = len(B)
	if not nA or not nB:
		return Enrichment(a,b,nA,nB,0,0,1)
	A.set_chromsizes("hg19")
	obs = len(A.intersect(B, u=True))
	dist = monte_carlo(B, A)
	exp = numpy.mean(dist)
	p_value = len([x for x in dist if x > obs]) / float(len(dist))
	if p_value > 0.5:
		p_value = 1 - p_value
	return Enrichment(a, basename(b), nA, nB, obs, exp, p_value)

"""
def annotation(x,y):
	for result in x.intersect(y,wa=True,wb=True):
		yield [x.name,y.name] + list(result)
"""

