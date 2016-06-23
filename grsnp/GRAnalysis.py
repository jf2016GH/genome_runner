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

class GRAnalysis:

	def __init__(self, fois_path, gfs_path, bg_path, outdir, job_name="", root_data_dir="", organism="",id="default",
				 console_output=False, cur_progress = None, print_progress = False):

		# setup logging
		self.logger = logging.getLogger(__name__)
		outpath = os.path.join(outdir, 'gr_log.txt')
		formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
		fh = logging.FileHandler(outpath)
		fh.setFormatter(formatter)
		self.logger.addHandler(fh)
		self.logger.setLevel(logging.INFO)
		self.progress_outpath = os.path.join(outdir,".prog")

		# Read in the paths
		fois = [line for line in self.read_lines(fois_path) if not line.endswith(".tbi")]
		gfs = [line for line in self.read_lines(gfs_path) if not line.endswith(".tbi")]
		# check if there are spaces in invalid parts of the file name
		invalid_names = utils.validate_filenames(fois + gfs + [bg_path])
		if len(invalid_names) != 0:
			self.logger.error("The following file(s) have invalid file names:\n" + "\n".join(invalid_names))
			self._write_progress("ERROR: Files have invalid filenames. See log file. Terminating run. See Analysis Log.",
							0)
			raise IOError("The following file(s) have invalid file names:\n" + "\n".join(invalid_names))
		if bg_path.endswith(".tbi"):
			self.logger.error("Background has invalid extension (.tbi). Terminating run.")
			self._write_progress("ERROR: Background has invalid extension (.tbi). Terminating run. See Analysis Log.",
							0)
			raise IOError("The following file(s) have invalid file names:\n" + "\n".join(invalid_names))

		# pre-process the FOIs
		fois = utils.preprocess_fois(fois, outdir, root_data_dir, organism)
		if len(fois) == 0:
			self.logger.error('No valid FOIs to supplied')
			self._write_progress("ERROR: No valid FOI files supplied. Terminating run. See Analysis Log.", 0)
			raise IOError("'No valid FOIs to supplied'")

		# pre-process the GFs and the background
		bg_path = utils.preprocess_gf_files([bg_path], root_data_dir, organism)[0]
		gfs = utils.preprocess_gf_files(gfs, root_data_dir, organism)

		# initialize class variables
		self.fois = fois
		self.gfs = gfs
		self.bg_path = bg_path
		self.outdir = outdir
		self.job_name = job_name
		self.root_data_dir = root_data_dir
		self.organism = organism
		self.id = id
		self.console_output = console_output
		self.cur_prog = 0 # int of the current progress
		self.max_prog = 0
		self.print_prog = print_progress


	def _write_progress(self, line, progress):
		"""Saves the current progress to the progress file
		progress: is a Project object which contains outpath, current value, and maximum value
		"""
		if not self.cur_prog:
			return
		if self.progress_outpath:
			dict_progress = {"status": line, "curprog": self.cur_prog, "progmax": self.max_prog}
			with open(progress.outpath, "wb") as progfile:
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


class GREnrichment(GRAnalysis):
	def __init__(self,fois, gfs, bg_path, outdir, job_name="", root_data_dir="", organism=""):
		GRAnalysis.__init__(self, fois, gfs, bg_path, outdir, job_name, root_data_dir, organism)
		self.detailed_outpath = os.path.join(outdir,"detailed.txt")
		self.matrix_outpath =  os.path.join(outdir, "matrix_PVAL.txt")
		self.matrix_sor_outpath = os.path.join(outdir, "matrix_OR.txt")

	def run_chisquare(self):
		pass

	def run_binomial(self):
		self._run_enrichment("binomial")

	def run_montecarlo(self):
		self._run_enrichment("chisquare")


	def _run_enrichment(self,stat_test):
		track_descriptions = []

		decriptions_path = os.path.join(self.root_data_dir, "grsnp_db", self.organism, "gf_descriptions.txt")
		if os.path.exists(decriptions_path):
			track_descriptions = [x.split("\t") for x in open(decriptions_path).read().split("\n") if x != ""]

		self.logger.info("Enrichment analysis started")
		# Validate FOIs against background. Also get the size of the background (n_bgs)
		foi_bg, good_fois = self.check_background_foi_overlap(self.bg_path, self.fois, progress=self.cur_prog)
		self.write_output("\t".join(map(base_name, good_fois)) + "\n", self.matrix_outpath)
		self.write_output("\t".join(map(base_name, good_fois)) + "\n", self.matrix_sor_outpath)
		self.write_output("\t".join(
			['foi_name', 'foi_obs', 'n_fois', 'bg_obs', 'n_bgs', 'odds_ratio', "ci_lower", "ci_upper",
			 'shrunken_odds_ratio', str(stat_test) + '_ p_val']) + "\n", self.detailed_outpath)
		self.cur_prog.current, self.cur_prog.max = 0, len(self.gfs)

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
				"###" + base_name(gf) + "\t" + self.get_score_strand_settings(gf) + "\t" + self._get_description(base_name(gf),
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

	def get_bgobs(self,bg, gf, root_data_dir, organism, progress=None):
		''' Check if pre-calculated GF and background overlap data exist.
		If they do not, it manually calculates them.
		'''
		# get the grsnp_db_[filt] folder
		filt_grsnp_db = gf.replace(root_data_dir, "").lstrip("/").split("/")[0]
		bkg_overlap_path = os.path.join(root_data_dir, filt_grsnp_db, organism, 'bkg_overlaps.gr')

		# See if pre-calculated values exist
		if os.path.exists(bkg_overlap_path):
			data = open(bkg_overlap_path).read().split("\n")
			data = [x.split("\t") for x in data if x != ""]
			d_gf = [x[1] for x in data if os.path.join(root_data_dir, x[0]) == gf and x[1] != ""]

			if len(d_gf) != 0:
				bg_obs = [x.split(":")[1] for x in d_gf[0].split(",") if x.split(":")[0] == os.path.basename(bg)]
				if len(bg_obs) != 0:
					self.logger.info("Pre-calculated values found for background and {} ".format(base_name(gf)))
					return bg_obs[0]
		# manually get overlap values
		self.logger.info("Calculating overlap stats on background and {}".format(base_name(gf)))
		if progress:
			self._write_progress("Calculating overlap stats on background and {}".format(base_name(gf)), progress)
		result = self.get_overlap_statistics(gf, [bg])
		try:
			result = int(result[0]["intersectregions"])
		except Exception, e:
			result = None
			self.logger.error(traceback.format_exc())
		return result

	def _output_p_value(self,foi_obs, n_fois, bg_obs, n_bgs, foi_path, gf_path, background_path, detailed_outpath,
					   stat_test=None):
		"""Return the shrunken odds-ratio and signed p-value of all FOIs against the GF so they can be written to
		matrix files. Outputs stats to the detailed results file.
		"""
		foi_name = base_name(foi_path)
		gf_name = base_name(gf_path)
		sign, pval, odds_ratio, shrunken_or, ci_lower, ci_upper = self._calculate_p_value_odds_ratio(foi_obs, n_fois, bg_obs,
																							   n_bgs, foi_name, gf_path,
																							   stat_test=stat_test,
																							   background_path=background_path,
																							   run_files_dir=self.outdir,
																							   progress=progress)
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
										 strpval])) + "\n", detailed_outpath)

		if pval < 1E-307:
			# set to value obtained from sys.float_info.min_10_exp
			pval = 1E-306
		return [sign * pval, shrunken_or]

