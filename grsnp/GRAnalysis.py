import collections
import gzip
import string
import os
import logging
import subprocess
import grsnp.analysis_util as utils
from grsnp.analysis_util import base_name as base_name
import traceback
import grsnp.dbcreator_util as db_utils
import json
import sys
import scipy
import scipy.stats
import math
import numpy as np

class GRAnalysis:

	def __init__(self, fois_path, gfs_path, bg_path, outdir, job_name="", root_data_dir="", organism="",job_id="default",
				 console_output=False, print_progress = False):
		try:
			# setup logging
			self.logger = logging.getLogger(__name__)
			outpath = os.path.join(outdir, 'gr_log.txt')
			formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
			fh = logging.FileHandler(outpath)
			fh.setFormatter(formatter)
			self.logger.handlers = [] # remove filehandler for previous runs if they exist
			self.logger.addHandler(fh)
			self.logger.setLevel(logging.INFO)
			self.progress_outpath = os.path.join(outdir,".prog")

			self.outdir = outdir
			self.job_name = job_name
			self.root_data_dir = root_data_dir
			self.organism = organism
			self.job_id = job_id
			self.console_output = console_output
			self.cur_prog = 0 # int of the current progress
			self.max_prog = 0
			self.print_prog = print_progress
			# Read in the paths
			fois = [line for line in utils.read_lines(fois_path) if not line.endswith(".tbi")]
			gfs = [line for line in utils.read_lines(gfs_path) if not line.endswith(".tbi")]
			# check if there are spaces in invalid parts of the file name
			invalid_names = self.validate_filenames(fois + gfs + [bg_path])
			if len(invalid_names) != 0:
				self._write_progress("ERROR: Files have invalid filenames. See log file. Terminating run. See Analysis Log.",
								0)
				raise IOError("Invalid file names")
			if bg_path.endswith(".tbi"):
				self.logger.error("Background has invalid extension (.tbi). Terminating run.")
				self._write_progress("ERROR: Background has invalid extension (.tbi). Terminating run. See Analysis Log.",
								0)
				raise IOError("Invalid file names")

			# pre-process the GFs and the background
			bg_path = self.preprocess_gf_files([bg_path])[0]
			gfs = self.preprocess_gf_files(gfs)

			# pre-process the FOIs
			fois = self.preprocess_fois(fois)
			if len(fois) == 0:
				self.logger.error('No valid FOIs to supplied')
				self._write_progress("ERROR: No valid FOI files supplied. Terminating run. See Analysis Log.", 0)
				raise IOError("'No valid FOIs to supplied'")

			# initialize class variables
			self.fois = fois
			self.gfs = gfs
			self.bg_path = bg_path


		except Exception, e:
			self.logger.error(traceback.print_exc())
			raise e
	def _write_progress(self, line, progress):
		"""Saves the current progress to the progress file
		progress: is a Project object which contains outpath, current value, and maximum value
		"""
		if not self.cur_prog:
			return
		if self.progress_outpath:
			dict_progress = {"status": line, "curprog": self.cur_prog, "progmax": self.max_prog}
			with open(self.progress_outpath, "wb") as progfile:
				progfile.write(json.dumps(dict_progress))
		if self.print_prog:
			print line

	def get_progress(self):
		"""returns the progress from the progress file
		"""
		path = os.path.join(self.outdir, ".prog")
		if os.path.exists(path):
			return open(path).read()
		else:
			return ""

	# Writes the output to the file specified.  Also prints to console if console_output is set to true
	def write_output(self, content, outpath=None):
		if outpath:
			with open(outpath, "a") as out:
				out.write(content)
		if self.console_output:
			print >> sys.stderr, content

	def _get_description(self, gf, track_descriptions):
		''' Get the GF feature description
		'''
		desc = [x[1] for x in track_descriptions if x[0] == gf and len(x[0]) > 1]
		if len(desc) is not 0:
			return desc[0]
		else:
			return "No Description"

	def validate_filenames(self, file_paths):
		''' Checks if there are spaces before or after the file extension.
		EX. 'dir1/dir2/test .bed' is not valid. 'dir1/dir2/test.bed is valid.
		'dir1/dir2/test.bed .gz' is not valid.
		'''
		invalid = []
		for file in file_paths:
			for t in os.path.basename(file).split("."):
				if len(t.strip()) != len(t):
					# there are spaces before or after the '.'. Add the file to the list of invalids.
					self.logger.error("Cannot have space before/in file extension: {}".format(file))
					invalid.append(os.path.basename(file))
		files_wo_ext = [base_name(x) for x in file_paths]
		file_duplicates = [x for x, y in collections.Counter(files_wo_ext).items() if y > 1]
		for f in file_duplicates:
			self.logger.error("{} exists multiple times. (i.e. 'foi.txt.gz' has same basename as 'foi.txt')".format(f))
			invalid.append(f)
		return invalid

	def preprocess_fois(self, fois):
		processed_fois = []
		output_dir = os.path.join(self.outdir, 'processed_fois')
		# Sort the fois
		out = ""
		data_dirs = [os.path.join(self.root_data_dir, 'grsnp_db', self.organism),
					 os.path.join(self.root_data_dir, 'custom_data')]
		try:
			for f in fois:
				if len([x for x in data_dirs if x in f]) == 0:
					# extract file if it is gzipped.
					unzipped_f = f
					if f.endswith('.gz'):
						out = subprocess.Popen(["gunzip {}".format(f)], shell=True, stdout=subprocess.PIPE,
											   stderr=subprocess.PIPE)
						out.wait()
						# filepath without the .gz extension
						unzipped_f = ".".join(unzipped_f.split(".")[:-1])
					# copy the FOI to the output_dir
					out_fname = os.path.split(unzipped_f)[1]
					if not out_fname.endswith('.bed'):
						out_fname = utils.base_name(out_fname) + ".bed"
					out_f = os.path.join(output_dir, out_fname)  # the processed foi file
					if not os.path.exists(output_dir):
						os.makedirs(output_dir)
					out = subprocess.Popen(['cp {} {}'.format(unzipped_f, out_f)], shell=True, stdout=subprocess.PIPE,
										   stderr=subprocess.PIPE)
					out.wait()
					# remove the header from the files
					db_utils.remove_headers(out_f)

					# Check if items are rsIDs. Convert to bed coordinates if they are
					if utils.validate_rsids(out_f):
						# check if a file exists in the database for rsID conversion and construct the path to it
						rsid_path = os.path.join(self.root_data_dir, 'custom_data', 'rsid_conversion', self.organism)
						if not os.path.exists(rsid_path):
							self.logger.error(
								'rsID conversion not available for this self.organism. Feature set {} removed'.format(f))
							continue
						files = [x for x in os.listdir(rsid_path) if
								 os.path.isfile(os.path.join(rsid_path, x)) and x.endswith('.bed')]
						# if conversion files found, perform conversion
						if len(files) > 0:
							rsid_path = os.path.join(rsid_path, files[0])
							script = """join {} {} -1 1 -2 4 -o 2.1 -o 2.2 -o 2.3 -o 2.4 -o 2.5 -o 2.6 | sed 's/\ /\t/g' > {}.temp""".format(
								out_f, rsid_path, out_f)
							out = subprocess.Popen([script], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
							out.wait()
							tmp_er = out.stderr.read()
							if tmp_er != "":
								self.logger.error(tmp_er)
								raise Exception(tmp_er)
							# we remove the original out_f FOI file and replace with the out_f.temp created with the join command above
							os.remove(out_f)
							out = subprocess.Popen(['cp {} {}'.format(out_f + '.temp', out_f)], shell=True,
												   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
							out.wait()
							os.remove(out_f + '.temp')
						else:
							self.logger.error('rsID conversion not available for this self.organism. Analysis terminated.')
							return []

						# sort the FOI and bgzip
					db_utils.sort_convert_to_bgzip(out_f, out_f + ".gz")
					script = "tabix -f " + out_f + ".gz"
					out = subprocess.Popen([script], shell=True, stdout=subprocess.PIPE)
					out.wait()
					# check if processed FOI file is empty
					if os.stat(out_f + ".gz")[6] != 0:
						processed_fois.append(out_f + ".gz")
					else:
						self.logger.error("{} is empty. Removing.")
				else:
					processed_fois.append(f)
				# print "{} is part of the database".format(f)
		except Exception, e:
			self.logger.error("Error while processing the FOIs")
			self.logger.error(traceback.format_exc())
			raise Exception(e)
		return processed_fois


	def preprocess_gf_files(self, files_paths):
		''' Used to preprocess the GF files and the background file.
		Returns gzipped file paths
		'''
		processed_files = []
		data_dirs = [os.path.join(self.root_data_dir, 'grsnp_db', self.organism),
					 os.path.join(self.root_data_dir, 'custom_data')]
		try:
			for out_f in files_paths:
				# if the file is not uploaded by the user, do not preprocess
				if len([x for x in data_dirs if x in out_f]) == 0:
					# unzip the file
					unzipped_f = out_f
					if out_f.endswith('.gz'):
						out = subprocess.Popen(["gunzip {}".format(out_f)], shell=True, stdout=subprocess.PIPE,
											   stderr=subprocess.PIPE)
						out.wait()
						# filepath without the .gz extension
						unzipped_f = ".".join(unzipped_f.split(".")[:-1])

					out_fname = os.path.split(unzipped_f)[1]
					if not out_fname.endswith('.bed'):
						out_fname = utils.base_name(out_fname) + ".bed"
					# replace the ".txt" extension with ".bed"
					output_dir = os.path.split(unzipped_f)[0]
					out_bed_f = os.path.join(output_dir, out_fname)  # the processed file
					db_utils.remove_headers(unzipped_f)
					# perform rsid conversion
					if utils.validate_rsids(unzipped_f):
						# check if a file exists in the database for rsID conversion and construct the path to it
						rsid_path = os.path.join(self.root_data_dir, 'custom_data', 'rsid_conversion', self.organism)
						if not os.path.exists(rsid_path):
							self.logger.error('rsID conversion not available for this self.organism. Feature set {} removed'.format(
								unzipped_f))
							continue
						files = [x for x in os.listdir(rsid_path) if
								 os.path.isfile(os.path.join(rsid_path, x)) and x.endswith('.bed')]
						# if conversion files found, perform conversion
						if len(files) > 0:
							rsid_path = os.path.join(rsid_path, files[0])
							# sort the RSID
							script = """sort -k1,1 -o {} {}""".format(unzipped_f, unzipped_f)
							out = subprocess.Popen([script], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
							out.wait()
							# join the RSID with the SNP data in custom_data
							script = """join {} {} -1 1 -2 4 -o 2.1 -o 2.2 -o 2.3 -o 2.4 -o 2.5 -o 2.6 | sed 's/\ /\t/g' > {}.temp""".format(
								unzipped_f, rsid_path, unzipped_f)
							out = subprocess.Popen([script], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
							out.wait()
							tmp_er = out.stderr.read()
							if tmp_er != "":
								self.logger.error(tmp_er)
								raise Exception(tmp_er)
							# we remove the original unzipped_f FOI file and replace with the unzipped_f.temp created with the join command above
							out = subprocess.Popen(['cp {} {}'.format(unzipped_f + '.temp', out_bed_f)], shell=True,
												   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
							out.wait()
							os.remove(unzipped_f + '.temp')
						else:
							self.logger.error('rsID conversion not available for this self.organism. Analysis terminated.')
							return []

					# sort the FOI
					db_utils.sort_convert_to_bgzip(out_bed_f, out_bed_f + ".gz")
					# check if processed FOI file is empty
					if os.stat(out_bed_f + ".gz")[6] != 0:
						# add the processed file (the ".bed" file) to the list.
						processed_files.append(out_bed_f + ".gz")
					else:
						self.logger.error("{} is empty. Removing.")

				else:
					processed_files.append(out_f)
					# print "{} is part of the database".format(out_f)
		except Exception, e:
			self.logger.error("Error while processing the GF/background files")
			self.logger.error(traceback.format_exc())
			raise Exception(e)
		return processed_files

	def _chunks(l, n):
		n = max(1, n)
		return [l[i:i + n] for i in range(0, len(l), n)]


class GRAnnotation(GRAnalysis):
	def __init__(self,fois, gfs, bg_path, outdir, job_name="", root_data_dir="", organism="",job_id="default",
				 console_output = False,print_progress = False):
		GRAnalysis.__init__(self,fois, gfs, bg_path, outdir, job_name=job_name, root_data_dir=root_data_dir,
							organism=organism,job_id=job_id, console_output = console_output, pring_progress=print_progress)

	def run_annotation(self):
		try:
			self.logger.info("Annotation started{}".format(self.job_id))
			annot_outdir = os.path.join(self.outdir, "annotations")
			if not os.path.exists(annot_outdir): os.mkdir(annot_outdir)
			self.cur_prog, self.max_prog = 0, len(self.fois)
			for f in self.fois:
				self._write_progress("Running Annotation Analysis for {}.".format(base_name(f)), self.cur_prog)
				self.logger.info("Running annotation analysis for {}".format(base_name(f)))
				for i, g in enumerate(utils._chunks(self.gfs, 100)):
					with open(os.path.join(annot_outdir, base_name(f) + str(i) + ".txt"), "wb") as wr:
						anot = self._get_annotation(f, g).split("\n")
						anot[0] = anot[0].replace("Region\t\t", "Region\t")
						wr.write("Region" + "\t" + "\t".join(base_name(x) for x in reversed(anot[0].split("\t")[
																							1:])) + "\tTotal")  # annotationAnalysis column order is reverse of input order
						for ind, a in enumerate(anot[1:]):
							if a.strip() != "":
								cur_row = a.split("\t")
								wr.write("\n" + str(ind) + "|" + "\t".join(
									cur_row + [str(sum([int(x) for x in cur_row[1:] if x != ""]))]))
				self.cur_prog += 1
			self.cur_prog, self.max_prog = 1, 1
			self._write_progress("Anotation finished", self.cur_prog)
			self.logger.info("Annotation finished for {}".format(self.job_id))
		except Exception, e:
			self.logger.error(traceback.format_exc())
			raise e

	def _get_annotation(self,foi, gfs):
		"""
		fois: list of FOI filepath
		gfs: filepaths for GF
		"""
		# use temporary files instead of piping out to console because large amounts of output to console can cause deadlock
		# this creates unique random file names
		tmp_path = utils.get_tmp_file('grsnptmp')
		tmp_error_path = utils.get_tmp_file('grsnperrortmp')

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
			if tmp_er != "": self.logger.error(tmp_er)
			if tmp[:6] == "ERROR:":
				self.logger.error(tmp[7:])
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
			self.logger.error(traceback.format_exc())
			raise e
		return tmp


class GREnrichment(GRAnalysis):
	def __init__(self,fois, gfs, bg_path, outdir, job_name="", root_data_dir="", organism="",job_id="default",console_output = False,print_progress = False):
		GRAnalysis.__init__(self, fois, gfs, bg_path, outdir, job_name, root_data_dir, organism,job_id=job_id,console_output=console_output,print_progress=print_progress)
		self.detailed_outpath = os.path.join(outdir,"detailed.txt")
		self.matrix_outpath =  os.path.join(outdir, "matrix_PVAL.txt")
		self.matrix_sor_outpath = os.path.join(outdir, "matrix_OR.txt")

	def run_chisquare(self):
		try:
			self._run_enrichment("chisquare")
		except Exception, e:
			self.logger.error(traceback.format_exc())
			raise e

	def run_binomial(self):
		try:
			self._run_enrichment("binomial")
		except Exception, e:
			self.logger.error(traceback.format_exc())
			raise e

	def run_montecarlo(self,num_mc):
		try:
			self._run_enrichment("montecarlo_" + str(num_mc))
		except Exception, e:
			self.logger.error(traceback.format_exc())
			raise e

	def _run_enrichment(self,stat_test):
		track_descriptions = []

		decriptions_path = os.path.join(self.root_data_dir, "grsnp_db", self.organism, "gf_descriptions.txt")
		if os.path.exists(decriptions_path):
			track_descriptions = [x.split("\t") for x in open(decriptions_path).read().split("\n") if x != ""]

		self.logger.info("Enrichment analysis started for {}".format(self.job_id))
		# Validate FOIs against background. Also get the size of the background (n_bgs)
		foi_bg, good_fois = self._check_background_foi_overlap(self.bg_path, self.fois, progress=self.cur_prog)
		self.write_output("\t".join(map(base_name, good_fois)) + "\n", self.matrix_outpath)
		self.write_output("\t".join(map(base_name, good_fois)) + "\n", self.matrix_sor_outpath)
		self.write_output("\t".join(
			['foi_name', 'foi_obs', 'n_fois', 'bg_obs', 'n_bgs', 'odds_ratio', "ci_lower", "ci_upper",
			 'shrunken_odds_ratio', str(stat_test) + '_ p_val']) + "\n", self.detailed_outpath)
		self.cur_prog, self.max_prog = 0, len(self.gfs)

		# check if any good fois exist after background filtering
		if len(good_fois) == 0:
			self.logger.error('No valid FOIs to supplied')
			self._write_progress("ERROR: No valid FOI files supplied. Terminating run. See Analysis Log.", self.cur_prog)
			return
		# remove old detailed enrichment result files if they exit
		enr_path = os.path.join(self.outdir, "enrichment")
		for f in good_fois:
			f_path = os.path.join(enr_path, base_name(f) + '.txt')
			if os.path.exists(f_path): os.remove(f_path)
		self._write_progress("Performing calculations on the background.", self.cur_prog)
		for gf in self.gfs:
			self._write_progress("Performing {} analysis for {}".format(stat_test, base_name(gf)), self.cur_prog)
			self.write_output(
				"###" + base_name(gf) + "\t" + self._get_score_strand_settings(gf) + "\t" + self._get_description(base_name(gf),
																												  track_descriptions) + "###" + "\n",
				self.detailed_outpath)
			res = self.get_overlap_statistics(gf, good_fois)

			# calculate bg_obs
			bg_obs = self.get_bgobs(gf)
			if bg_obs == None:
				self.logger.error("Skipping {}".format(gf))
				continue

			n_bgs = foi_bg[0]["indexregions"]

			# calculate the pvalues and output the matrix line for the current gf
			pvals, sors = [], []  # sors = shrunken odds-ratios
			for i in range(len(good_fois)):
				[pvalue, shrunken_or] = self._output_p_value(res[i]["intersectregions"], res[i]["queryregions"], bg_obs,
													   n_bgs, good_fois[i], gf, stat_test=stat_test)
				pvals.append(str(pvalue))
				sors.append(str(shrunken_or))

			# output the matrices file lines
			self.write_output("\t".join([base_name(gf)] + pvals) + "\n", self.matrix_outpath)
			self.write_output("\t".join([base_name(gf)] + sors) + "\n", self.matrix_sor_outpath)

			self.cur_prog += 1

		self.cur_prog, self.max_prog = 1, 1
		self._write_progress("Analysis Completed", self.cur_prog)
		self.logger.info("Analysis Completed {}".format(self.job_id))

	def _get_score_strand_settings(self, gf_path):
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

	def get_overlap_statistics(self,gf, fois):
		"""Returns a dictionary with indicating how many hits exist for each foi against the gf
		gf: filepath for GF
		fois: list of FOI filepaths
		"""
		results = []
		out = ""
		# use temporary files instead of piping out to console because large amounts of output to console can cause deadlock
		# this creates unique random file names
		tmp_path = utils.get_tmp_file('grsnptmp')
		tmp_error_path = utils.get_tmp_file('grsnperrortmp')
		tmp_file = open(tmp_path, 'wb')
		tmp_error_file = open(tmp_error_path, 'wb')

		try:
			# Runs overlapStatistics with preprocessed background stats if they exist
			out = subprocess.Popen(["overlapStatistics"] + [gf] + fois, stdout=tmp_file, stderr=tmp_error_file)
			out.wait()
			tmp_file.close()
			tmp_error_file.close()
			tmp = open(tmp_path).read()
			tmp_er = open(tmp_error_path).read()
			if tmp_er != "": self.logger.error(tmp_er)
			if tmp[:6] == "ERROR:":
				self.logger.error(tmp[7:])
				raise Exception(tmp)

			for x in tmp.split("\n")[1:]:
				if x != "":
					tmp = x.split("\t")
					foi_name, n, hit_count = os.path.split(tmp[0])[-1], tmp[2], tmp[3]
					results.append({"queryfile": foi_name, "queryregions": int(n), "intersectregions": int(hit_count),
									"indexregions": int(tmp[1])})
			# remove the temporary output files
			if os.path.exists(tmp_path): os.remove(tmp_path)
			if os.path.exists(tmp_error_path): os.remove(tmp_error_path)
		except Exception, e:
			if not tmp_file.closed: tmp_file.close()
			if not tmp_error_file.closed: tmp_error_file.close()
			# remove the temporary output files
			if os.path.exists(tmp_path): os.remove(tmp_path)
			if os.path.exists(tmp_error_path): os.remove(tmp_error_path)
			self.logger.error(traceback.format_exc())
			raise e
		return results

	def get_bgobs(self, gf):
		''' Check if pre-calculated GF and background overlap data exist.
		If they do not, it manually calculates them.
		'''
		# get the grsnp_db_[filt] folder
		filt_grsnp_db = gf.replace(self.root_data_dir, "").lstrip("/").split("/")[0]
		bkg_overlap_path = os.path.join(self.root_data_dir, filt_grsnp_db, self.organism, 'bkg_overlaps.gr')

		# See if pre-calculated values exist
		if os.path.exists(bkg_overlap_path):
			data = open(bkg_overlap_path).read().split("\n")
			data = [x.split("\t") for x in data if x != ""]
			d_gf = [x[1] for x in data if os.path.join(self.root_data_dir, x[0]) == gf and x[1] != ""]

			if len(d_gf) != 0:
				bg_obs = [x.split(":")[1] for x in d_gf[0].split(",") if x.split(":")[0] == os.path.basename(self.bg_path)]
				if len(bg_obs) != 0:
					self.logger.info("Pre-calculated values found for background and {} ".format(base_name(gf)))
					return bg_obs[0]
		# manually get overlap values
		self.logger.info("Calculating overlap stats on background and {}".format(base_name(gf)))
		self._write_progress("Calculating overlap stats on background and {}".format(base_name(gf)), self.cur_prog)
		result = self.get_overlap_statistics(gf, [self.bg_path])
		try:
			result = int(result[0]["intersectregions"])
		except Exception, e:
			result = None
			self.logger.error(traceback.format_exc())
		return result

	def _output_p_value(self,foi_obs, n_fois, bg_obs, n_bgs, foi_path, gf_path, stat_test=None):
		"""Return the shrunken odds-ratio and signed p-value of all FOIs against the GF so they can be written to
		matrix files. Outputs stats to the detailed results file.
		"""
		foi_name = base_name(foi_path)
		gf_name = base_name(gf_path)
		sign, pval, odds_ratio, shrunken_or, ci_lower, ci_upper = self._calculate_sor_pvall(foi_obs, n_fois, bg_obs,
																							n_bgs, foi_name, gf_path)
		if stat_test == "binomial":
			pval = self._calc_binom_pval(foi_obs, n_fois, bg_obs, n_bgs)
		elif stat_test.startswith("motecarlo"):
			num_mc = stat_test.split("_")[1]
			pval = self._calc_mc_pval(foi_obs, n_fois, gf_path, num_mc, odds_ratio)


		if sign == 1 or str(odds_ratio) == "inf":
			direction = "overrepresented"
		else:
			direction = "underrepresented"

		strpval = ""
		self.write_output("\t".join(map(str, [foi_name.rpartition('/')[-1], foi_obs, n_fois, bg_obs, n_bgs,
											  utils._format_num(odds_ratio),
											  utils._format_num(ci_lower),
											  utils._format_num(ci_upper),
											  utils._format_num(shrunken_or),
										 "%.2e" % pval if type(pval) != type("") else pval,
										 strpval])) + "\n", self.detailed_outpath)

		if pval < 1E-307:
			# set to value obtained from sys.float_info.min_10_exp
			pval = 1E-306
		return [sign * pval, shrunken_or]

	def _calculate_sor_pvall(self, foi_obs, n_fois, bg_obs, n_bgs, foi_name, gf_path, do_chi_square = True):
		"""Calculates the p-value,confidence intervals and the shrunken odds ratio.
		Returns [sign,pval,odds_ratio,shrunken_or,ci_lower,ci_upper]
		"""
		do_chi_square = True
		self._write_progress("Testing {}".format(foi_name), self.cur_prog)
		# Perform the chisquare test regardless of what stat_test is selected, we need the odds ratio
		bg_obs, n_bgs = int(bg_obs), int(n_bgs)
		ctable = [[foi_obs, n_fois - foi_obs],
				  [bg_obs - foi_obs, n_bgs - n_fois - (bg_obs - foi_obs)]]
		# Ensure there are no negative values in the ctable
		for i in ctable:
			for k in i:
				if k < 0:
					self.logger.warning(
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

		sor = self._calc_shrunke_or(odds_ratio, ctable)

		return [sor["sign"], pval, odds_ratio, sor["shrunken_or"], sor["ci_lower"], sor["ci_upper"]]

	def _calc_mc_pval(self, foi_obs, n_fois, gf_path, num_mc, odds_ratio):
		'''
		Calculates the pvalue using Monte Carlo (MC) simulations.
		:param foi_obs:
		:param n_fois:
		:param gf_path:
		:param num_mc: number of MC to run
		:param odds_ratio: Required to determine whether fois are under or overrepresentated ????
		:return: pvalue
		'''

		rndfoipath = os.path.join(self.outdir, 'mc.bed')
		# pow_mc states what starting power of 10 to check pvalue
		chunk_size, pow_mc, not_significant = 100, 2, False
		num_rnd_obs = []  # stores the number of rnd_snps that overlap for each mc

		# run the rnd_fois in groups against the GF (allows us to handle case of >10,000 MC simulations)
		for i_chunk in xrange(1, num_mc, chunk_size):
			if not_significant == True: break
			# only create the number of rnd_snps files needed (i.e for 14 mc with chunk of 10 we only want to create 4 files for last chunk)
			rnd_count = chunk_size if i_chunk + chunk_size < num_mc else num_mc - i_chunk + 1
			# Generate the random fois
			rnd_fois_paths = utils.generate_randomsnps(rndfoipath, self.bg_path, n_fois, rnd_count)
			#    _write_progress("Performing Monte Carlo {} of {}".format(i_chunk,num_mc), progress)
			# get overlap stats for random_features against the GF
			overlapstats = self.get_overlap_statistics(gf_path, rnd_fois_paths)
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
					self._write_progress(
						"Pval at {} runs calculated as {}. Numerator: {} Denominator: {}".format(i_chunk + i_res, pval,
																								 float(num_over) + 1,
																								 float(len(
																									 num_rnd_obs)) + 1),
						self.cur_prog)
					if pval >= 1 / float(pow(10, pow_mc)):
						# pval will never be significant stop doing Monte Carlo
						not_significant = True
						self._write_progress("Not significant. Stopping Monte carlo (P-value = {}".format(pval), self.cur_prog)
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
		return pval

	def _calc_binom_pval(self, foi_obs, n_fois, bg_obs, n_bgs):
		'''
		Calculate the pvalue usig the binomial distribution
		:param foi_obs:
		:param n_fois:
		:param bg_obs:
		:param n_bgs:
		:return: pvalue
		'''
		pval = scipy.stats.binom_test(foi_obs, n_fois, float(bg_obs) / n_bgs)
		return pval

	def _calc_shrunke_or(self, odds_ratio,ctable):
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
		sign = -1 if odds_ratio < 1 else 1
		return {"sign": sign, "shrunken_or": shrunken_or, "ci_lower": ci_lower, "ci_upper": ci_upper}

	def _check_background_foi_overlap(self, bg, fois, progress=None):
		""" Calculates the overlap of the FOIs with the background.
		Removes FOIs that are poorly formed with the background.
		"""
		if progress:
			self._write_progress("Validating FOIs against background", progress)
		good_fois = []
		if len(fois) == 0:
			return [[], []]
		# Runs overlapStatistics on background and FOIs
		foi_bg_stats = self.get_overlap_statistics(bg, fois)
		for f in foi_bg_stats:
			isgood = True
			foi_name, n_bgs, n_fois, foi_in = f["queryfile"], f["indexregions"], f["queryregions"], f[
				"intersectregions"]
			if n_fois < 5:
				isgood = False
				self.logger.warning("Number of SNPs in {} < 5. Removing it from analysis.".format(foi_name))
			elif n_bgs < n_fois:
				isgood = False
				self.logger.warning("Number of SNPs in {} > than in background. Removing it from analysis.".format(foi_name))
			if isgood:
				# ensure that overlapStatistics output filename with extension for queryFile field
				good_fois.append([x for x in fois if os.path.basename(x) == f["queryfile"]][0])
			if foi_in < n_fois:
				self.logger.warning(
					"{} out of {} {} SNPs are not a part of the background. P-value are unreliable. Please, include all SNPs in the background and re-run analysis.".format(
						n_fois - foi_in, n_fois, foi_name))
		return [foi_bg_stats, good_fois]
