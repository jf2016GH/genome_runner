#!/usr/bin/env python2
"""
Analysis of several FOI files against several GFs using Fisher's exact test. Best used for SNP set analysis, using whole SNP database as a spot background.
sys.argv[1] - text file with FOI file names. Include full, or relative path, if needed.
sys.argv[2] - text file with GF file names. Include full, or relative path, if needed.
sys.argv[3] - spot background file
"""

import argparse
import collections
import gzip
import math
import sys
import logging
from logging import FileHandler,StreamHandler
from bx.intervals.intersection import IntervalTree
from scipy.stats import hypergeom
import numpy
import scipy
import pdb
import os
import json

Interval = collections.namedtuple("Interval", "chrom,start,end")

# Logging configuration
logger = logging.getLogger()
logger = logging.getLogger('genomerunner.query')
hdlr = logging.FileHandler('genomerunner_server.log')
hdlr_std = StreamHandler()
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
# This line outputs logging info to the console
logger.addHandler(hdlr_std)
logger.setLevel(logging.INFO)


matrix_outpath = None
detailed_outpath = None
progress_outpath = None
console_output = False

def read_intervals(path, background=None, snp_only=False):
    with (gzip.open(path) if path.endswith(".gz") else open(path)) as h:
        intervals = []
        for i,line in enumerate(h):
            fields = line.strip().split("\t")
#            print fields
            chrom, start, end = fields[:3]
            start = int(start)
            end = int(end)
            if snp_only and ((end - start) > 1):
                logger.warning("\t"+" ".join([chrom,str(start),str(end)]) + " in {} is not a SNP. Shortening feature to be SNP".format(path))
                end = start + 1
            interval = Interval(chrom, int(start), int(end))
            if background and not background.query(interval):
                logger.warning("\t"+" ".join([chrom,str(start),str(end)]) + " in {} does not overlap with the background. Continuing...".format(path))
            intervals.append(interval)
    return intervals

class IntervalSet(object):
    def __init__(self, intervals, snp_only=False):
        self._trees = {}
        self._intervals = list(intervals)
        for i,iv in enumerate(intervals):
            if not iv.chrom in self._trees:
                self._trees[iv.chrom] = IntervalTree()
            self._trees[iv.chrom].add(iv.start, iv.end, i)

    def query(self, interval):
        t = self._trees.get(interval.chrom)
        if t:
            return t.find(interval.start, interval.end)

    def n_overlaps(self, intervals):
        return sum([1 if self.query(iv) else 0 for iv in intervals])

    def __len__(self):
        return len(self._intervals)

def p_value(gf, fois, bgs, foi_name):
    "Return the signed log10 p-value of intersection with the given interval set."
    foi_obs = gf.n_overlaps(fois) # number of FOIs overlapping with a GF
    bg_obs = gf.n_overlaps(bgs) # number of spot bkg overlapping with a GF
    n_fois, n_bgs = len(fois), len(bgs)
    ctable = [[foi_obs, n_fois-foi_obs],
              [bg_obs-foi_obs,n_bgs-(bg_obs-foi_obs)]]
    odds_ratio, pval = scipy.stats.fisher_exact(ctable)
    sign = 1 if (odds_ratio < 1) else -1
    write_output("\t".join(map(str, [foi_name.rpartition('/')[-1], foi_obs, n_fois, bg_obs, n_bgs, "%.2f" % odds_ratio, "%.2f" % pval])) + "\n",detailed_outpath)
    try:        
        return sign * math.log10(pval)        
    except ValueError as e:
        logger.error(e)
        return sign * sys.float_info.min_10_exp


def run_hypergeom(fois, gfs, bg_path,outpath):    
    global detailed_outpath,matrix_outpath, progress_outpath, curprog, progmax
    detailed_outpath =  os.path.join(outpath, "detailed.gr")
    matrix_outpath = os.path.join(outpath,"matrix.gr")
    progress_outpath = os.path.join(outpath,".prog")
    fois = read_lines(fois)
    gfs = read_lines(gfs)

    bg = read_intervals(bg_path)
    bg_iset = IntervalSet(bg)
    foi_sets = dict((path,read_intervals(path, snp_only=True, background=bg_iset)) for path in fois)

    write_output("\t".join(fois)+"\n", matrix_outpath)
    write_output("\t".join(['foi_name', 'foi_obs', 'n_fois', 'bg_obs', 'n_bgs', 'odds_ratio', 'p_val']) + "\n",detailed_outpath)
    curprog,progmax = 0,len(gfs)    
    for gf in gfs:
        curprog += 1
        _write_progress("Performing Hypergeometric analysis for {}".format(gf))
        gf_iset = IntervalSet(read_intervals(gf))
        write_output(gf+"\n",detailed_outpath)
        [str(p_value(gf_iset, foi_sets[foi], bg, foi)) for foi in fois]
        write_output("\t".join([gf] + [str(p_value(gf_iset, foi_sets[foi], bg, foi)) for foi in fois])+"\n",matrix_outpath)
    _write_progress("Analysis Completed")

# Writes the output to the file specified.  Also prints to console if console_output is set to true
def write_output(content,outpath=None):
    if outpath:
        with open(outpath,"a") as out:
            out.write(content)
    if console_output:
        print >> sys.stderr, content

# Collect lines from a file into a list
def read_lines(path):
    elems = []
    with open(path) as h:
        for line in h:
            if line.strip():
                elems.append(line.strip())
    return elems

def _write_progress(line):
    """Saves the current progress to the progress file
    """
    global progress_outpath
    if progress_outpath:
        global curprog, progmax
        progress = {"status": line, "curprog": curprog,"progmax": progmax}
        with open(progress_outpath,"wb") as progfile:
            print "writePROGRESS", progress
            progfile.write(json.dumps(progress))

if __name__ == "__main__":
    console_output = True
    parser = argparse.ArgumentParser(description="Create a matrix of hypergeometric p-values for genomic intersections.")
    parser.add_argument("fois", nargs=1, help="Text file with FOI file names (SNPs only).") 
    parser.add_argument("gfs", nargs=1, help="Text file with GF file, gzipped.") 
    parser.add_argument("bg_path", nargs=1, help="Path to spot background file (SNPs only).")
    args = parser.parse_args()

    fois = read_lines(args.fois[0])
    gfs = read_lines(args.gfs[0])

    bg = read_intervals(args.bg_path[0])
    bg_iset = IntervalSet(bg)
    foi_sets = dict((path,read_intervals(path, snp_only=True, background=bg_iset)) for path in fois)

    write_output("\t".join(fois)+"\n", matrix_outpath)
    write_output("\t".join(['foi_name', 'foi_obs', 'n_fois', 'bg_obs', 'n_bgs', 'odds_ratio', 'p_val'])+"\n",detailed_outpath)
    for gf in gfs:
        gf_iset = IntervalSet(read_intervals(gf))
        write_output(gf+"\n",detailed_outpath) 
        write_output("\t".join([gf] + [str(p_value(gf_iset, foi_sets[foi], bg, foi)) for foi in fois])+"\n",matrix_outpath)
