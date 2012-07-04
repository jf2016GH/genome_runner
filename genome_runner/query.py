#!/usr/bin/env python
"""
GenomeRunner "lite" edition.
"""
import sys, operator, os, numpy
from pybedtools import BedTool
from collections import namedtuple

Enrichment = namedtuple("Enrichment",
		["A","B","nA","nB","observed","expected","p_value"])
def enrichment_category(self):
	return "nonsig"
Enrichment.category = enrichment_category

def monte_carlo(ref, query, n=10):
	result = []
	for i in range(n):
		shuffled = query.shuffle(g="data/hg19.genome").sort()
		overlaps = shuffled.intersect(ref,u=True)
		result.append(len(overlaps))
	return result

def make_filter(name, score, strand):
	def filter(interval):
		if name and name != interval.name:
			return False
		if score and score > interval.score:
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
	A = BedTool(a).filter(flt)
	B = BedTool(b).filter(flt)
	A.set_chromsizes("hg19")
	obs = len(A.intersect(B))
	dist = monte_carlo(B, A)
	exp = numpy.mean(dist)
	p_value = len([x for x in dist if x > obs]) / float(len(dist))
	if p_value > 0.5:
		p_value = 1 - p_value
	return Enrichment(a, b, len(A), len(B),
			obs, exp, p_value)

"""
def annotation(x,y):
	for result in x.intersect(y,wa=True,wb=True):
		yield [x.name,y.name] + list(result)
"""

