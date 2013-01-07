#!/usr/bin/env python
import sys, operator, os, numpy
import datetime
from numpy import random as rand
from pybedtools import BedTool
from pybedtools.featurefuncs import midpoint as pb_midpoint
import pybedtools
from collections import namedtuple
import cPickle
import logging
from logging import FileHandler,StreamHandler
from os.path import basename
import json
import copy
import traceback  as trace
from scipy.stats import uniform,kstest,hypergeom
import pprint,argparse,pdb,scipy, math


logger = logging.getLogger('genomerunner.query')
hdlr = logging.FileHandler('genomerunner_server.log')
hdlr_std = StreamHandler()
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
# This line outputs logging info to the console
logger.addHandler(hdlr_std)
logger.setLevel(logging.INFO)
res_path = "results"
DEBUG = True


# This class represents an Enrichment analysis result
# Lists of Enrichment objects are serialized to a Python Pickle
# file when an analysis is complete	
_Enrichment = namedtuple("Enrichment",
		["A","B","nA","nB","observed","expected","p_value","obsprox","expprox","pybed_p_value",
		"pybed_expected","jaccard_observed","jaccard_p_value",
		"jaccard_expected","proximity_p_value","kolmogorov_p_value","hypergeometric_p_value"])
_Enrichment_Par = namedtuple("Enrichment_Par","a b A B Background n flt genome genome_fn organism obs background")
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


# TODO implement remove_invalid and report results to us 10 	er.  
def enrichment(id,a, b,background, organism,name=None, score=None, strand=None, n=10, run=[]):
	"""Perform enrichment analysis between two BED files.

	a - path to Feature of Interest BED file
	b - path to Genomic Feature BED file
	n - number of Monte-Carlo iterations
	"""
	write_debug("START",True)
	r = {}
	e = _Enrichment_Par
	e.a,e.b,e.organism,e.n,e.Background, e.background = a,b,organism,n,None,background

	if not background is None:
		e.Background = BedTool(background)

	e.A = BedTool(str(e.a))
	e.B = BedTool(str(e.b))
	e.genome = pybedtools.get_chromsizes_from_ucsc(e.organism)
	e.genome_fn = pybedtools.chromsizes_to_file(e.genome)

	e.organism = str(e.organism)
	flt = make_filter(name,score,strand)
	e.B.filter(flt).saveas()
	e.nA = len(e.A)
	e.nB = len(e.B)
	if not e.nA or not e.nB:
		return Enrichment(e.a,basename(e.b),e.nA,e.nB,0,0,1,0,0,1,0,0,1,0)

	e.A.set_chromsizes(e.genome)
	e.B.set_chromsizes(e.genome)
	e.obs = len(e.A.intersect(e.B, u=True))
	# This is the Monte-Carlo step.  If custom background present, it is used
	if 'pvalue' in run:
		logger.info("Running Monte Carlo ({}): (id={})".format(b,id))
		write_progress(id, "Running Monte Carlo {}".format(b))
		r.update(run_montecarlo(_Enrichment_Par))
	else:
		r['p_value'], r['exp'] = "NA","NA"
		logger.info("Skipping Monte Carlo ({}): (id={})".format(b,id))
	## Uncomment to print global parameters in debug file
	#write_debug("Global parameters",a = e.a,b=e.b,A= e.A,B=e.B,background = e.background,
	#			n=e.n,flt = e.flt,genome = e.genome,genome_fn = e.genome_fn,organism = e.organism,obs = e.obs)
	
	# expected caluclated using pybed method CANNOT use custom background
	if 'pybedtool' in run:
		logger.info("Running Random Intersections ({}): (id={})".format(b,id))
		write_progress(id, "Running Random Intersections: {0}".format(b))
		r.update( run_pybedtool(_Enrichment_Par))
	else:
		r['pybedp_value'], r['pybed_exp'] = "NA","NA"
		logger.info("Skipping Random Intersections")

	# epected calculated using jaccard method
	if 'jaccard' in run:
		logger.info("Running Jaccard ({}): (id={})".format(e.b,id))
		write_progress(id, "Running Jaccard {}".format(e.b))
		r.update( run_jaccard(_Enrichment_Par))
	else:
		r['jaccardp_value'], r['jaccard_obs'],r['jaccard_exp'] = "NA","NA","NA"
		logger.info("Skipping Jaccard ({}): (id={})".format(b,id))
	
	# run kolmogorov-smornov test
	if 'kolmogorov' in run:
		logger.info("Running Kolmogorov-Smornov {} (id={})".format(b,id))
		write_progress(id, "Running Kolmogorov-Smornov{}".format(b))
		r.update( run_kolmogorov(_Enrichment_Par))
	else:
		r['kol_smor_p_value'] = "NA"
		logger.info("Skipping Kolmogorov-Smornov {} (id={})".format(b,id))
	
	# run proximity analysis
	if 'proximity' in run:
		logger.info("Running proximity {} (id={})".format(b,id))
		write_progress(id, "Running proximity analysis{}".format(b))
		r.update( run_proximity(_Enrichment_Par))
	else:
		logger.info( "Skipping Proximity")
		r['obsprox'],r['expprox'],r['proximityp_value']="NA","NA", "NA" 

	# run hypergeometric distrubtion analysis
	if 'hypergeometric' in run:
		write_progress(id,"Running")
		logger.info("Running hypergeometric analysis {} (id={})".format(b,id))
		r.update(run_hypergeometric(_Enrichment_Par))
	else:
		logger.info("Skipping hypergeometric")
		r['hypergeometric_p_value'] = "NA"
	## Uncomment to print global parameters in debug file
	#write_debug("Global parameters",a = e.a,b=e.b,A= e.A,B=e.B,background = e.background,
	#		n=e.n,flt = e.flt,genome = e.genome,genome_fn = e.genome_fn,organism = e.organism,obs = e.obs)
	# the order of these arguments IS IMPORTANT
	return Enrichment(e.a, basename(e.b), e.nA, e.nB, e.obs, r['exp'], r['p_value'],r['obsprox'],\
			r['expprox'],r['pybedp_value'],r['pybed_exp'],r['jaccard_obs'],r['jaccardp_value'],\
			r['jaccard_exp'],r['proximityp_value'],r['kol_smor_p_value'],r['hypergeometric_p_value'])


def run_montecarlo(Enrichment_Par):
	e = Enrichment_Par
	if e.background is not None:
		dist = [len(e.A.shuffle(genome=e.organism,chrom=True,incl=e.background).intersect(e.B, u=True)) for i in range(e.n)]
	else:
		dist = [len(e.A.shuffle(genome=e.organism,chrom=True).intersect(e.B, u=True)) for i in range(e.n)]

	exp = numpy.mean(dist)
	# gave p_value a value here so that it doesn't go out of scope, is this needed?
	p_value = 'NA'
	if exp == e.obs or (exp == 0 and e.obs == 0):
		p_value =1
	else:
		p_value = len([x for x in dist if x > e.obs]) / float(len(dist))
		p_value = min(p_value, 1 - p_value)		

	write_debug(run_montecarlo.__name__,False,Enrichment_Par= str(e),p_value= p_value, observed = e.obs, expected =exp)
	return {"exp": exp,"p_value": p_value}

def run_pybedtool(Enrichment_Par):
	e = Enrichment_Par
	pybeddist = e.A.randomintersection(e.B,iterations=e.n,shuffle_kwargs={'chrom': True})
	pybeddist = list(pybeddist)
	pybed_exp = numpy.mean(pybeddist)
	pybedp_value = 'NA'
	if pybed_exp == e.obs or (pybed_exp == 0 and e.obs == 0):
		pybedp_value =1
	else:
		pybedp_value = len([x for x in pybeddist if x > e.obs]) / float(len(pybeddist))
		pybedp_value = min(pybedp_value,1-pybedp_value)
	return { "pybedp_value":pybedp_value, "pybed_exp": pybed_exp} 

def run_jaccard(Enrichment_Par):
	e = Enrichment_Par
	# cuts out chrom,chromStart, and chromEnd
	A2 = e.A.cut([0,1,2])
	B2 = e.B.cut([0,1,2])	   
	genome_fn = pybedtools.chromsizes_to_file(e.genome)
	resjaccard = A2.random_jaccard(B2,genome_fn=genome_fn,iterations=e.n,
			shuffle_kwargs={'chrom':True})
	jaccard_dist = resjaccard[1]
	jaccard_obs = resjaccard[0]
	jaccard_exp = numpy.mean(resjaccard[1])
	jaccardp_value = 'NA'
	if jaccard_exp == jaccard_obs or (jaccard_exp  == 0 and jaccard_obs == 0):
		jaccardp_value =1
	else:
		jaccardp_value = len([x for x in jaccard_dist if x > jaccard_obs]) / float(len(jaccard_dist))
		jaccardp_value = min(jaccardp_value,1-jaccardp_value)
	return {"jaccardp_value": jaccardp_value, "jaccard_obs": jaccard_obs, "jaccard_exp": jaccard_exp} 

def run_kolmogorov(Enrichment_Par):
	e = Enrichment_Par
	return {"kol_smor_p_value": genome_tri_corr(e.A,e.B,e.genome)}


def run_proximity(Enrichment_Par):
	#stores the means of the distances for the MC
	e = Enrichment_Par 
	expall =[]
	for i in range(e.n):
		tmp = shuffle(e.A,e.Background,e.organism).closest(e.B,d=True)
		# get the distances
		for t in tmp:
			expall.append(t[-1])
	# calculate the overal expected distance
	expall.append(numpy.mean(numpy.array(expall,float)))
	# calculate the expected mean for all of the runs
	expprox = numpy.mean(numpy.array(expall,float))	
	# proximety analysis for observed
	tmp = e.A.closest(e.B,d=True)
	obsall = []
	for t in tmp:
		obsall.append(t[-1])
	obsprox = numpy.mean(numpy.array(obsall,float))
	proximityp_value = len([x for x in expall if x > obsprox]) / float(len(expall))
	proximityp_value = min(proximityp_value,1-proximityp_value)
	return {"expprox": expprox, "proximityp_value": proximityp_value,"obsprox":obsprox}


def run_hypergeometric(Enrichment_Par):
	"""  Runs the hypergeometric test.  Requires a background. Returns 'NA' if no background
	is provided.
	"""
	e = Enrichment_Par
	gf,foi,bg,hypergeomp_value = e.A,e.B,e.Background,"NA"
	if e.Background is not None:
		foi_obs = len(foi.intersect(gf, u = True)) # number of FOIs overlapping with a GF
		bg_obs = len(bg.intersect(gf, u = True)) # number of spot bkg overlapping with a GF
		rnd_obs = (len(foi)*bg_obs/len(bg)) # Mean of hypergeometric distiribution
		# pdb.set_trace()
		if foi_obs == rnd_obs: # No difference
			hypergeomp_value =  0
		elif foi_obs < rnd_obs: # Underrepresentation		
			pval = scipy.stats.fisher_exact([[foi_obs,bg_obs],[len(foi)-foi_obs,len(bg)-bg_obs]],alternative='less')[1]
			hypergeomp_value = math.log10(pval) # Calculating cdf hypergeometric distribution, not adding "-" to signify underrepresentation
		elif foi_obs > rnd_obs:	# Overrepresentation		
			pval = scipy.stats.fisher_exact([[foi_obs,bg_obs],[len(foi)-foi_obs,len(bg)-bg_obs]],alternative='greater')[1]
			hypergeomp_value = -math.log10(pval) # Calculating sf hypergeometric distribution
		else: # No difference
			hypergeomp_value = 0
			# print "e.obs : %d, rand.obs :			%d, len(back) : %d, back.obs : %d, len(e.A) : %d, hypergeomp_value : %r" % (eobs, randobs, len(back), backobs, len(eA), hypergeomp_value)
	return {"hypergeometric_p_value": hypergeomp_value}

def shuffle(bedtool, background,organism):
	""" Accepts a file path to a bed file and a background file
	Random intervals are generate from the intervals in the background
	file.  An attempt is made to make the new random intervals the same size
	as the original intervals. One new interval is generated per interval in the
	bed file.

	If no background is provided, The pybedtool internal shuffle function is called
	and the entire genome is used as background
	"""
	A, B =  bedtool,background
	print "test "
	if B is not None:
		rand_a = ""
		for a in A:
			a_len = a.length
			r = rand.random_integers(0,len(B)-1)
			# if cur A length is greater than random background inteval
			# the new randA is set to be same size as background inteval
			if  a_len > B[r].length:
				rand_a += "\t".join(B[r][:4]) +"\t" + "\t".join(a[4:])+"\n"
			else:
				randstart = rand.random_integers(B[r].start,B[r].end - a_len)
				rand_a += "\t".join([B[r].chrom,str(randstart),str(randstart+a_len),
					B[r].name,
					"\t".join(a[4:])]) + "\n"	
		rand_A = BedTool(rand_a,from_string=True)
		return rand_A
	else:
		return A.shuffle(genome=organism,chrom=True)

def genome_tri_corr(A_bedtool, B_bedtool,genome):
	''' Based on Genometric Correlation algorithm at http://genometricorr.sourceforge.net/
		Genome file is a dictionary of chroms and chrom lengths.  Can be retrieved from ucsc
		by calling pybedtool.get_chromsizes_from_ucsc('hg19'). Where hg19 is the 
		assembly name.
		Is strand specific.
	'''
	# generate the midpoint intervals
	midpoints_A = A_bedtool.each(_midpoint).saveas()
	midpoints_B = B_bedtool.each(_midpoint).saveas()
	midpoints_B = midpoints_B.sort().saveas()
	midpoints_B = _fill_ends(midpoints_B,genome).saveas()

	closest_A = midpoints_A.closest(midpoints_B,d=True,s=True)
	d_numor_A = [x[-1] for x in closest_A]


	# gets the closest midpoint upstream
	closest_up = midpoints_A.closest(midpoints_B,D="b",io=True,iu=True,s=True)
	# gets the closest midpoint downstream
	closest_down = midpoints_A.closest(midpoints_B,D="b",io=True,id=True,s=True)


	d_denom_A = []
	#calculates the denominator values for each d_i
	map(lambda u,d: d_denom_A.append(abs(int(u[-1])) + abs(int(d[-1]))), closest_up,closest_down)

	D = []
	# calculates the D for each foi in A 
	map(lambda n,d: D.append(float(n)/float(d)), d_numor_A,d_denom_A)
	if len(D) != 0:
		# performs the kolmogorov-smornov test
		rv=uniform(low=0, high=1, scale=0.5)
		results = kstest(D,rv.cdf)
		return results[1] 
	else:
		return "NA"
			



def _fill_ends(B_bedtool,genome):
	''' Creates midpoints for the chrom start and chrom end of the genomic feature.'''

	# save a new copy of the bed tool
	B = B_bedtool.saveas()
	# get all the chromosome
	chroms = set()
	for i in range(0,len(B)):
		chroms.add(B[i].chrom)
	B_extremes = ""
	# create a midpoint for the start and end of the chromosome for each strand
	for c in chroms:
		# gets the strands that contain gf on the current chrom
		strands = set()
		map(lambda x: strands.append(x.strand),B.filter(lambda x: x.chrom == c))
		if len(strands) == 0: strands.add(".")

		for s in strands:
			c_start = genome[c][0]
			c_end = genome[c][1]
			B_extremes += "{}\t{}\t{}\t{}\t\t\n".format(c,c_start,c_start+1,s)
			B_extremes += "{}\t{}\t{}\t{}\t\t\n".format(c,c_end,c_end+1,s)
	B_out = pybedtools.BedTool(B_extremes,from_string=True)
	return B.cat(B_out,post_merge=False).saveas()



def _midpoint(interval):
	''' Uses pybedtools midpoint function to calculate the midpoint of the interval
	Sets the endpoint to start + 1 to conformt to the bed format standard.
	'''
	

	#print "before midpoint: " + str(interval)
	i = pb_midpoint(interval)
	#print "after  midpoint: " + str(i)
	i.end = i.start + 1
	#print "after endpoint +1" + str(i)
	return i

	

def generate_background(foipath,gfpath,background):
	"""accepts a background filepath
	generate a background and returns as a pybedtool.
	Replaces the chrom fields of the foi and the gf with the interval
	id from the background.
	"""
	bckg = background
	bckgnamed = "" 
	interval = 0 

	#inserts a unique interval id into the backgrounds name field
	for b in bckg:
		bckgnamed +=  "\t".join(b[:3])+'\t{}\t'.format(interval) + "\t".join(b[4:]) + "\n"
		interval += 1
	bckg = BedTool(bckgnamed,from_string=True)
	foi = BedTool(str(foipath))
	gf = BedTool(str(gfpath))
	# get the interval names from the background that the gf intersects with
	gf = bckg.intersect(gf)
	gfnamed = ""

	# insert the interval id into the chrom field of the gf and creates a new bedtool
	for g in gf:
		gfnamed += '{}\t'.format(g.name) + "\t".join(g[1:]) + "\n"
		#print "GFNAMED: " + str(g)
	gf = BedTool(gfnamed,from_string=True)
	#print "GFBEDTOOL: " + str(g)

	# inserts the interval id into the chrom column of the foi and creates a new bedtool
	foi = bckg.intersect(foi)
	foinamed = ""
	for f in foi:
		foinamed += '{}\t'.format(f.name) + "\t".join(f[1:])+"\n" 
		#print "FOINAMED: " + str(f)
	foi = BedTool(foinamed,from_string=True)

	#print "FOIBEDTOOL: " + str(f)
	bckgnamed = ""
	for b in bckg:
		bckgnamed += '{}\t'.format(b.name) + "\t".join(b[1:])+"\n"
	bckg = BedTool(bckgnamed,from_string=True)
	# converts the background to a genome dictionary
	chrstartend = [(g.start,g.end) for g in bckg]
	background = dict(zip([g.chrom for g in bckg],chrstartend))
	return {"foi": foi,"gf":gf,"background":background}


	

	run_pvalue=False,run_pybedtool=False,run_jaccard=False,run_proximity=False,run_kolmogorov=False
def run_enrichments(id, f, gfeatures,background, niter, name, score, strand,organism,run):
	"""
	Run one FOI file (f) against multiple GFs, then 
	save the result to the "results" directory. Returns list of enrichment results.
	"""
	# sets up logging for the run
	if not os.path.exists(res_path): os.makedirs(res_path)
	hdlr_id_file = logging.FileHandler(os.path.join(res_path,str(id)+".log"))
	logger.addHandler(hdlr_id_file)
	try:
		enrichments = []
		# these are progress values that are written to the progress file
		global curprog 
		curprog = 0
		global progmax
		progmax = len(gfeatures)		
		header = True
		for gf in gfeatures:
			write_progress(id, "RUNNING ENRICHMENT ANALYSIS FOR: {}".format(gf))
			e = enrichment(id,f,gf,background,organism,name,score,strand,niter,run)
			enrichments.append(e)
			path = os.path.join(res_path, str(id))
			with open(path, "w") as strm:
				cPickle.dump(enrichments, strm)
			curprog += 1
			write_results_astext(id,e,header)
			header = False
		write_progress(id, "FINISHED")
	except Exception, e:
		logger.error(e)
		logger.error(trace.format_exc())
		write_progress(id,"ERROR: The run crashed: {}".format(e))
	return enrichments



def write_progress(id,line):
	"""Saves the current progress to the progress file
	"""
	global curprog
	global progmax
	path = os.path.join(res_path,str(id)+".prog")
	progress = {"status": line, "curprog": curprog,"progmax": progmax}
	with open(path,"wb") as progfile:
		progfile.write(json.dumps(progress))

# enrichment_results is the namedtuple generated from the enrichment
def write_results_astext(id,enrichment_results,header=False):
	"""Writes the results of the enrichment analysis into a plain text file.
	"""
	fields = {"A": "FOI_Name","B": "GF_Name","nA": "#FOI","nB": "#GF","observed": "observed",
		"expected":"expected","p_value":"p_value","obsprox":"obsprox","expprox":"expprox",
		"pybed_p_value":"pybed_p_value","pybed_expected":"pybed_expected","jaccard_observed":"jaccard_observed"
		,"jaccard_p_value":"jaccard_p_value","jaccard_expected":"jaccard_expected","proximity_p_value":"proximity_p_value"
		,"kolmogorov_p_value":"kolmogorov_p_value","hypergeometric_p_value":"hypergeometric_p_value"}
	
	outpath =  os.path.join(res_path,str(id)) + ".txt"
	
	if header:	
		with open(outpath,"w") as w:
			line = ""		
			for f in _Enrichment._fields:
				line += fields[f]+"\t" 
			line = line[:-1]
			w.write(line+"\n")
	with open(outpath,"a") as w:
		line = ""
		for i in range(len(enrichment_results)):
			line += str(enrichment_results[i]) + "\t"
		line = line[:-1]
		w.write(line + "\n")

		
def get_progress(id):
	"""returns the progress from the progress file
	"""

	path = os.path.join(res_path,str(id) + ".prog")
	if os.path.exists(path):
		return open(path).read()
	else:
		return ""


def write_debug(fun_name,header = False,**kwargs):
	if DEBUG:
		with open("Debug.log","a") as w:
			if header:
				curtime = "\nDEBUG\t" + str(datetime.datetime.now().strftime("%y-%m-%d %H:%M")) + "\n"
				w.write(curtime)
			else:
				w.write("Variable values for\t %s \n" % fun_name)
				for key, value in kwargs.iteritems():
					w.write("%s =\t %s\n" % (key,value))


# id, f, gfeatures,background, niter, name, score, strand,organism,run
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Runs Enrichment analysis in GenomeRunner')
	parser.add_argument('--jobid','-i',
		help='The file name to output the results to',default="")
	parser.add_argument('--foifile','-a', 
		help='File containing a list of features of interest file paths. Required')
	parser.add_argument('--gffile','-b', 
		help='File containing a list of genomic feature file paths. Required')	
	parser.add_argument('--background','-k', 
		help='File containging background in bed format.  Uses organism genome as background by default',default = None)
	parser.add_argument('--organism','-g', 
		help='The UCSC code of the organism to be downloaded (example: hg19 (human))',default = 'hg19')
	parser.add_argument('--mcnum','-n', 
		help='The number of Monte Carlo runs to perform',default = 10)
	parser.add_argument('--run','-r', dest="run", action='append', 
		help='Test to run.  Repeat to run multiple test. Available tests: pvalue,pybedtool,jaccard,kolmogorov,proximity,hypergeometric')
	parser.add_argument('--name','-e', help='Name to filter by.',default = None)
	parser.add_argument('--strand','-s', help='Strand to filter by',default = None)
	parser.add_argument('--score','-c', help='Score to filter by',default = None)


	args = vars(parser.parse_args())

	res_path = "../results_cmd"

	if not os.path.exists(res_path):
		os.makedirs(res_path)
	if args['foifile'] is None or args['gffile'] is None:
		print "foifile and gffile parameters may not be blank.  Use --help for details."
		sys.exit()
	gfs = open(args['gffile'],'r').readlines()
	gfs = [line.strip() for line in gfs]
	header = True
	for foi in open(args['foifile'],'r').readlines():
		foi = foi.strip()
		result_id = args['jobid']+ "_" + os.path.splitext(os.path.basename(foi))[0]

		score = args['score']
		if not score is None:
			score = int(args['score']) 

		results = run_enrichments(result_id,foi,gfs,args['background'],int(args['mcnum']),args['name'],score,args['strand'],
						args['organism'],args['run'])
		header = False
