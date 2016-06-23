#!/usr/bin/env python2
from __future__ import division

import argparse
import json
import logging
import math
import numpy as np
import os
import scipy
import subprocess
import sys
import traceback
from scipy.stats import hypergeom

from grsnp.analysis_util import base_name, get_tmp_file, _zip_run_files, \
	validate_filenames, generate_randomsnps


run_files_dir  = ""


# This line outputs logging info to the console
console_output = False
print_progress = False


def calculate_p_value_odds_ratio(foi_obs, n_fois, bg_obs, n_bgs, foi_name, gf_path, stat_test=None,
								background_path=None, progress=None, run_files_dir=None):
	"""Calculates the p-value,confidence intervals and the shrunken odds ratio.
	Returns [sign,pval,odds_ratio,shrunken_or,ci_lower,ci_upper]
	"""
	_write_progress("Testing {}".format(foi_name), progress)
	# Perform the chisquare test regardless of what stat_test is selected, we need the odds ratio
	bg_obs, n_bgs = int(bg_obs), int(n_bgs)
	ctable = [[foi_obs, n_fois - foi_obs],
			  [bg_obs - foi_obs, n_bgs - n_fois - (bg_obs - foi_obs)]]
	# Ensure there are no negative values in the ctable
	do_chi_square = True
	for i in ctable:
		for k in i:
			if k < 0:
				logger.warning(
					"Cannot calculate p-value for {} and {}. Is the background too small? foi_obs {}, n_fois {}, bg_obs {}, n_bgs {}".format(
						base_name(gf_path), foi_name, foi_obs, n_fois, bg_obs, n_bgs))
				return [1, 1, 1, 1, 1, 1]
	# # ??? if sample too small, then perform fisher exact test
	#        if k < 5:
	#            do_chi_square = False
	# check for zeros and add 0.5 if one of the cells is 0
	if ctable[0][0] == 0 or ctable[0][1] == 0 or ctable[1][0] == 0 or ctable[1][1] == 0:
		ctable[0][0] += 0.5
		ctable[0][1] += 0.5
		ctable[1][0] += 0.5
		ctable[1][1] += 0.5

	if do_chi_square:
		chi_result = scipy.stats.chi2_contingency(ctable)
		pval = chi_result[1]
		odds_ratio = float(ctable[0][0] * ctable[1][1]) / (ctable[0][1] * ctable[1][0])
	else:
		odds_ratio, pval = scipy.stats.fisher_exact(ctable)
	# Adjustments of outliers
	if odds_ratio == 0.0:
		odds_ratio = sys.float_info.min
	if np.isinf(odds_ratio):
		odds_ratio = sys.float_info.max
	# # If p-value is insignificant, so is odds ratio
	#    if pval == 1.0:
	#        odds_ratio = 1

	# calculate the shrunken odds ratio
	log_or = scipy.log(odds_ratio)
	conf_coe = 1.96  # the confidence coefficient of a standard norm dist
	# calculate the standard error
	se = math.sqrt(1.0 / ctable[0][0] + 1.0 / ctable[1][0] + 1.0 / ctable[0][1] + 1.0 / ctable[1][1])
	# calculate the upper and lower confidence interval
	ci_upper = scipy.exp(log_or + conf_coe * se)
	ci_lower = scipy.exp(log_or - conf_coe * se)
	# Precaution against CI overflow
	if np.isinf(ci_upper):
		ci_upper = sys.float_info.max
	if ci_lower == 0.0:
		ci_lower = sys.float_info.min
	# shrunken_or is the ci (either upper or lower) that is closest to
	if odds_ratio < 1:
		ci_array = [odds_ratio, ci_upper if ci_upper < 1 else odds_ratio]
		ci_index = scipy.array(ci_array).argmax()
		shrunken_or = ci_array[ci_index]
	elif odds_ratio > 1:
		ci_array = [ci_lower if ci_lower > 1 else odds_ratio, odds_ratio]
		ci_index = scipy.array(ci_array).argmin()
		shrunken_or = ci_array[ci_index]
	else:
		shrunken_or = 1
	# if ci_lower<1 and ci_upper>1:
	#     shrunken_or,odds_ratio = 1,1
	# else:
	# find which value is closer to 1
	# ci_index = scipy.array([[abs(math.log(ci_lower)),abs(math.log(ci_upper))]]).argmin()
	# shrunken_or = [ci_lower,ci_upper][ci_index]

	## If a different stat_test is selected, perform that test now, and replace the p-value
	## note we will still use the odds ratio calculated by the chi-square test
	if stat_test == "binomial":
		pval = scipy.stats.binom_test(foi_obs, n_fois, float(bg_obs) / n_bgs)

	# monte carlo is passed as 'montecarlo_[number_of_simulations]'
	elif stat_test.startswith("montecarlo"):
		num_mc = int(stat_test.split("_")[1])
		rndfoipath = os.path.join(run_files_dir, 'mc.bed')
		# pow_mc states what starting power of 10 to check pvalue
		chunk_size, pow_mc, not_significant = 100, 2, False
		num_rnd_obs = []  # stores the number of rnd_snps that overlap for each mc

		# run the rnd_fois in groups against the GF (allows us to handle case of >10,000 MC simulations)
		for i_chunk in xrange(1, num_mc, chunk_size):
			if not_significant == True: break
			# only create the number of rnd_snps files needed (i.e for 14 mc with chunk of 10 we only want to create 4 files for last chunk)
			rnd_count = chunk_size if i_chunk + chunk_size < num_mc else num_mc - i_chunk + 1
			# Generate the random fois
			rnd_fois_paths = generate_randomsnps(rndfoipath, background_path, n_fois, rnd_count)
			#    _write_progress("Performing Monte Carlo {} of {}".format(i_chunk,num_mc), progress)
			# get overlap stats for random_features against the GF
			overlapstats = get_overlap_statistics(gf_path, rnd_fois_paths)
			for i_res, res in enumerate(overlapstats):
				if not_significant == True: break
				# get the rnd_obs
				num_rnd_obs.append(float(res["intersectregions"]))
				# check if we are at 10^(pow_mc)th result
				if i_chunk + i_res == pow(10, pow_mc):
					# Count how many random snp sets have more observed than foi_obs
					num_over = sum([1 for rnd_i in num_rnd_obs if rnd_i >= foi_obs])
					# calculate pvalue
					pval = (float(num_over) + 1) / (float(len(num_rnd_obs)) + 1)
					# Calculate depletion p-values, if necessary
					if odds_ratio < 1:
						pval = 1 - pval
					# pval = (1/float(pow(10,pow_mc)) - sys.float_info.min) if pval == 0 else pval
					_write_progress(
						"Pval at {} runs calculated as {}. Numerator: {} Denominator: {}".format(i_chunk + i_res, pval,
																								 float(num_over) + 1,
																								 float(len(
																									 num_rnd_obs)) + 1),
						progress)
					if pval >= 1 / float(pow(10, pow_mc)):
						# pval will never be significant stop doing Monte Carlo
						not_significant = True
						_write_progress("Not significant. Stopping Monte carlo (P-value = {}".format(pval), progress)
					pow_mc += 1
			for f in rnd_fois_paths:
				os.remove(f)

				# Calculate Monte carlo p-value
		# Count how many random snp sets have more observed than foi_obs
		num_over = sum([1 for rnd_i in num_rnd_obs if rnd_i >= foi_obs])
		# calculate pvalue
		pval = (float(num_over) + 1) / (float(len(num_rnd_obs)) + 1)
		# Calculate depletion p-values, if necessary
		if odds_ratio < 1:
			pval = 1 - pval
		pval = 1 / float(pow(10, pow_mc)) if pval == 0 else pval

	sign = -1 if odds_ratio < 1 else 1
	return [sign, pval, odds_ratio, shrunken_or, ci_lower, ci_upper]


def _chunks(l, n):
	n = max(1, n)
	return [l[i:i + n] for i in range(0, len(l), n)]


def get_annotation(foi, gfs):
	"""
	fois: list of FOI filepath
	gfs: filepaths for GF
	"""
	# use temporary files instead of piping out to console because large amounts of output to console can cause deadlock
	# this creates unique random file names
	tmp_path = get_tmp_file('grsnptmp')
	tmp_error_path = get_tmp_file('grsnperrortmp')

	tmp_file = open(tmp_path, 'wb')
	tmp_error_file = open(tmp_error_path, 'wb')
	try:
		out = subprocess.Popen(["annotationAnalysis"] + [foi] + gfs, stdout=tmp_file,
							   stderr=tmp_error_file)  # TODO enable ["--print-region-name"]
		out.wait()
		tmp_file.close()
		tmp_error_file.close()
		tmp = open(tmp_path).read()
		tmp_er = open(tmp_error_path).read()
		if tmp_er != "": logger.error(tmp_er)
		if tmp[:6] == "ERROR:":
			logger.error(tmp[7:])
			raise Exception(tmp)
		# remove the temporary output files
		if os.path.exists(tmp_path): os.remove(tmp_path)
		if os.path.exists(tmp_error_path): os.remove(tmp_error_path)

	except Exception, e:
		if not tmp_file.closed: tmp_file.close()
		if not tmp_error_file.closed: tmp_error_file.close()
		# remove the temporary output files
		if os.path.exists(tmp_path): os.remove(tmp_path)
		if os.path.exists(tmp_error_path): os.remove(tmp_error_path)
		logger.error(traceback.format_exc())
		raise e
	return tmp






def _write_head(content, outpath):
	f = front_appender(outpath)
	f.write(content)
	f.close()


def check_background_foi_overlap(bg, fois, progress=None):
	""" Calculates the overlap of the FOIs with the background.
	Removes FOIs that are poorly formed with the background.
	"""
	if progress:
		_write_progress("Validating FOIs against background", progress)
	good_fois = []
	if len(fois) == 0:
		return [[], []]
	# Runs overlapStatistics on background and FOIs
	foi_bg_stats = get_overlap_statistics(bg, fois)
	for f in foi_bg_stats:
		isgood = True
		foi_name, n_bgs, n_fois, foi_in = f["queryfile"], f["indexregions"], f["queryregions"], f["intersectregions"]
		if n_fois < 5:
			isgood = False
			logger.warning("Number of SNPs in {} < 5. Removing it from analysis.".format(foi_name))
		elif n_bgs < n_fois:
			isgood = False
			logger.warning("Number of SNPs in {} > than in background. Removing it from analysis.".format(foi_name))
		if isgood:
			# ensure that overlapStatistics output filename with extension for queryFile field
			good_fois.append([x for x in fois if os.path.basename(x) == f["queryfile"]][0])
		if foi_in < n_fois:
			logger.warning(
				"{} out of {} {} SNPs are not a part of the background. P-value are unreliable. Please, include all SNPs in the background and re-run analysis.".format(
					n_fois - foi_in, n_fois, foi_name))
	return [foi_bg_stats, good_fois]



def get_score_strand_settings(gf_path):
	''' Parses the gf_path and determines if gf is filtered by score and/or strand.
	'''
	str_strand, str_scorethresh = "Strand: Both", "Score threshold: NA"
	gfsplit = gf_path.split("/grsnp_db_")
	if len(gfsplit) == 2:
		str_score_strand = gfsplit[-1].split("/")[0].split("_")
		for s in str_score_strand:
			if s.isdigit():
				str_scorethresh = "Score threshold: " + s
			else:
				str_strand = "Strand: " + s
	return str_strand + "\t" + str_scorethresh


def run_hypergeom(fois, gfs, bg_path, outdir, job_name="", zip_run_files=False, bkg_overlaps_path="", root_data_dir="",
				  run_annotation=True, pct_score="", organism="", stat_test=None):
	valid_stat_tests = ["chisquare", "binomial"]
	if not os.path.exists(os.path.normpath(outdir)): os.mkdir(os.path.normpath(outdir))
	run_files_dir = outdir
	try:
		track_descriptions = []

		decriptions_path = os.path.join(root_data_dir, "grsnp_db", organism, "gf_descriptions.txt")
		if os.path.exists(decriptions_path):
			track_descriptions = [x.split("\t") for x in open(decriptions_path).read().split("\n") if x != ""]
		# set output settings

		if stat_test not in valid_stat_tests and not stat_test.startswith("montecarlo_"):
			logger.error("Valid p-value test not selected. Terminating run.")
			_write_progress("ERROR: Valid p-value test not selected. Terminating run. See Analysis Log.", curr_prog)
			print "ERROR: Valid p-value test not selected. Terminating run. See Analysis Log."
			return


		# Validate FOIs against background. Also get the size of the background (n_bgs)
		foi_bg, good_fois = check_background_foi_overlap(bg_path, fois, progress=curr_prog)
		write_output("\t".join(map(base_name, good_fois)) + "\n", matrix_outpath)
		write_output("\t".join(map(base_name, good_fois)) + "\n", matrix_sor_outpath)
		write_output("\t".join(
			['foi_name', 'foi_obs', 'n_fois', 'bg_obs', 'n_bgs', 'odds_ratio', "ci_lower", "ci_upper",
			 'shrunken_odds_ratio', str(stat_test) + '_ p_val',
			 'p_mod' if run_randomization_test else ""]) + "\n", detailed_outpath)
		curr_prog.current, curr_prog.max = 0, len(gfs)
		# check if any good fois exist after background filtering
		if len(good_fois) == 0:
			logger.error('No valid FOIs to supplied')
			_write_progress("ERROR: No valid FOI files supplied. Terminating run. See Analysis Log.", curr_prog)
			return
		# remove old detailed enrichment result files if they exit
		enr_path = os.path.join(run_files_dir, "enrichment")
		for f in good_fois:
			f_path = os.path.join(enr_path, base_name(f) + '.txt')
			if os.path.exists(f_path): os.remove(f_path)
		_write_progress("Performing calculations on the background.", curr_prog)
		for gf in gfs:
			current_gf = base_name(gf)
			_write_progress("Performing {} analysis for {}".format(stat_test, base_name(gf)), curr_prog)
			write_output(
				"###" + base_name(gf) + "\t" + get_score_strand_settings(gf) + "\t" + get_description(base_name(gf),
																									  track_descriptions) + "###" + "\n",
				detailed_outpath)
			res = get_overlap_statistics(gf, good_fois)

			# calculate bg_obs
			bg_obs = get_bgobs(bg_path, gf, root_data_dir, organism, progress=curr_prog)
			if bg_obs == None:
				logger.error("Skipping {}".format(gf))
				continue

			n_bgs = foi_bg[0]["indexregions"]

			# calculate the pvalues and output the matrix line for the current gf
			pvals, sors = [], []  # sors = shrunken odds-ratios
			for i in range(len(good_fois)):
				[pvalue, shrunken_or] = output_p_value(res[i]["intersectregions"], res[i]["queryregions"], bg_obs,
														n_bgs, good_fois[i], gf, bg_path, detailed_outpath,
														stat_test=stat_test, progress=curr_prog)
				pvals.append(str(pvalue))
				sors.append(str(shrunken_or))

			# output the matrices file lines
			write_output("\t".join([base_name(gf)] + pvals) + "\n", matrix_outpath)
			write_output("\t".join([base_name(gf)] + sors) + "\n", matrix_sor_outpath)

			curr_prog.current += 1
		if run_annotation:
			logger.info("Annotation started")
			annot_outdir = os.path.join(outdir, "annotations")
			if not os.path.exists(annot_outdir): os.mkdir(annot_outdir)
			curr_prog.current, curr_prog.max = 0, len(fois)
			for f in fois:
				_write_progress("Running Annotation Analysis for {}.".format(base_name(f)), curr_prog)
				logger.info("Running annotation analysis for {}".format(base_name(f)))
				for i, g in enumerate(_chunks(gfs, 100)):
					with open(os.path.join(annot_outdir, base_name(f) + str(i) + ".txt"), "wb") as wr:
						anot = get_annotation(f, g).split("\n")
						anot[0] = anot[0].replace("Region\t\t", "Region\t")
						wr.write("Region" + "\t" + "\t".join(base_name(x) for x in reversed(anot[0].split("\t")[
																							1:])) + "\tTotal")  # annotationAnalysis column order is reverse of input order
						for ind, a in enumerate(anot[1:]):
							if a.strip() != "":
								cur_row = a.split("\t")
								wr.write("\n" + str(ind) + "|" + "\t".join(
									cur_row + [str(sum([int(x) for x in cur_row[1:] if x != ""]))]))
				curr_prog.current += 1
			logger.info("Annotation finished")
		if zip_run_files:
			_write_progress("Preparing run files for download", curr_prog)
			_zip_run_files(outdir, job_name)
		curr_prog.current, curr_prog.max = 1, 1
		_write_progress("Analysis Completed", curr_prog)
		logger.info("Analysis Completed")
	except Exception, e:
		logger.error(traceback.print_exc())
		_write_progress("Run crashed. See end of log for details.", curr_prog)
		raise Exception(e)


class front_appender:
	'''
	Appends content to start of file.
	'''

	def __init__(self, fname, mode='a'):
		self.__write_queue = []
		self.__old_content = ""
		if mode == 'a':
			self.__old_content = open(fname).read()
		self.__f = open(fname, 'w')

	def write(self, s):
		self.__write_queue.append(s)

	def close(self):
		self.__f.writelines(self.__write_queue + [self.__old_content])
		self.__f.close()


def _load_minmax(path):
	data = {}
	if not os.path.exists(path):
		return data
	score = [x for x in open(path).read().split("\n") if x != ""]
	for s in score:
		name, min_max = s.split('\t')
		data[name] = min_max
	return data


def main():
	global detailed_outpath, progress_outpath, run_files_dir, console_output, print_progress
	print_progress = True
	parser = argparse.ArgumentParser(
		description="Enrichment analysis of several sets of SNPs (FOIs) files against several genomic features (GFs). Example: python analysis.py foi_full_names.txt gf_full_names.txt /path_to_background/snp137.bed.gz")
	parser.add_argument("fois", nargs=1, help="Text file with paths to FOI files (unless -p used). Required")
	parser.add_argument("gfs", nargs=1,
						help="Text file with pathrs to GF files (unless -p used). GF files may be gzipped. Required")
	parser.add_argument("bg_path", nargs=1, help="Path to background, or population of all SNPs. Required")
	parser.add_argument("--run_annotation", "-a", help="Run annotation analysis", action="store_true")
	parser.add_argument("--run_files_dir", "-r", nargs="?",
						help="Set the directory where the results should be saved. Use absolute path. Example: /home/username/run_files/.",
						default=os.getcwd())
	parser.add_argument("--pass_paths", "-p",
						help="Pass fois and gfs as comma separated paths. Paths are saved in .fois and .gfs file.",
						action="store_true")
	parser.add_argument("--data_dir", "-d", nargs="?", type=str,
						help="Set the directory containing the database. Required for rsID conversion. Use absolute path. Example: /home/username/db_#.##_#.##.####/.",
						default="")
	parser.add_argument('--organism', '-g', nargs="?",
						help="The UCSC code of the organism to use. Required for rsID conversion. Default: hg19 (human).",
						default="hg19")
	default_test = "chisquare"
	parser.add_argument('--stat_test', '-s', nargs="?",
						help="Select the statistical test to use for calculating P-values. Default: {}. Available: chisquare, binomial, montecarlo_[# of simulations]".format(
							default_test), default=default_test)
	args = vars(parser.parse_args())
	if args['organism'] is None:
		print "--organism cannot be blank"
		return None
	if args['run_files_dir'] is None:
		print "--run_files_dir cannot be blank"
		return None
	if args["pass_paths"]:
		gf = args["gfs"][0].split(",")
		foi = args["fois"][0].split(",")
		if not os.path.exists(args["run_files_dir"]) and args['run_files_dir'] != "":
			os.mkdir(args['run_files_dir'])
		# write out the passed gf and foi paths into .gfs and .fois files.
		args["gfs"][0], args["fois"][0] = os.path.join(args["run_files_dir"], ".gfs"), os.path.join(
			args["run_files_dir"], ".fois")
		with open(args["gfs"][0], 'wb') as writer:
			writer.write("\n".join(gf))
		with open(args["fois"][0], "wb") as writer:
			writer.write("\n".join(foi))
	run_hypergeom(args["fois"][0], args["gfs"][0], args["bg_path"][0], args["run_files_dir"], "", False, "",
				  args['data_dir'], args["run_annotation"], run_randomization_test=False, organism=args['organism'],
				  stat_test=args['stat_test'])


if __name__ == "__main__":
	main()


