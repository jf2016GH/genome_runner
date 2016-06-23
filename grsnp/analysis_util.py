import collections
import commands
import gzip
import os
import random
import string
import subprocess
import tarfile
import traceback
import zipfile

from grsnp import dbcreator_util as grsnp_util
from grsnp.analysis import logger
import datetime



def base_name(k):
	return os.path.basename(k).split(".")[0]

# Collect lines from a file into a list
def read_lines(path):
	elems = []
	with open(path) as h:
		for line in h:
			if line.strip():
				elems.append(line.strip())
	return elems

def _format_num(num):
	''' Sets format to be either scientific or float depending on num value
	'''
	if type(num) != type(""):
		if num > 100 or num < 0.01:
			return "%.2e" % num
		else:
			return "%.2f" % num
	else:
		return num


def get_tmp_file(prefix):
	tmp_path = prefix + "_" + ''.join(random.choice(string.lowercase + string.digits) for _ in range(32)) + '.tmp'
	while (os.path.exists(tmp_path)):
		tmp_path = prefix + "_" + ''.join(random.choice(string.lowercase + string.digits) for _ in range(32)) + '.tmp'
	return tmp_path


def _zip_run_files(outdir, id=""):
	'''
	File paths of FOIs and GFs as a list. Gathers all the files together in one zipped file
	'''

	# zip annotation result folder if it exists
	anno_dir = os.path.join(outdir, 'annotations')
	if os.path.exists(anno_dir):
		z_ano = zipfile.ZipFile(os.path.join(outdir, 'annotations.zip'), 'a')
		for f in os.listdir(anno_dir):
			z_ano.write(os.path.join(anno_dir, f), f)
		z_ano.close()

	# zip processed FOIs directory
	proc_dir = os.path.join(outdir, 'processed_fois')
	if os.path.exists(proc_dir):
		z_ano = zipfile.ZipFile(os.path.join(outdir, 'processed_fois.zip'), 'a')
		for f in os.listdir(proc_dir):
			z_ano.write(os.path.join(proc_dir, f), f)
		z_ano.close()

	f = open(os.path.join(outdir, "gr_log.txt"))
	f_log = f.read()
	f.close()
	path_settings, f_sett = os.path.join(outdir, ".settings"), ""
	if os.path.exists(path_settings):
		f = open(path_settings)
		f_sett = f.read() + "\n###LOG###\n"
		f.close()
	new_log_path = os.path.join(outdir, "gr_log.txt")
	new_log = open(new_log_path, 'wb')
	new_log.write(f_sett + f_log)
	new_log.close()
	tar_path = os.path.join(outdir, 'GR_{}.tar'.format(id))
	tar = tarfile.TarFile(tar_path, "a")
	output_files = [os.path.join(outdir, x) for x in os.listdir(outdir) if
					x.endswith(".txt") or x.endswith(".pdf") or x.endswith('.zip')]
	fls = output_files
	for f in fls:
		tar.add(f, os.path.basename(f))
	tar.close()
	tar_file = open(tar_path, 'rb')
	with gzip.open(tar_path + ".gz", "wb") as gz:
		gz.writelines(tar_file)
	tar_file.close()
	if os.path.exists(tar_path): os.remove(tar_path)


def validate_filenames(file_paths):
	''' Checks if there are spaces before or after the file extension.
	EX. 'dir1/dir2/test .bed' is not valid. 'dir1/dir2/test.bed is valid.
	'dir1/dir2/test.bed .gz' is not valid.
	'''
	invalid = []
	for file in file_paths:
		for t in os.path.basename(file).split("."):
			if len(t.strip()) != len(t):
				# there are spaces before or after the '.'. Add the file to the list of invalids.
				logger.error("Cannot have space before/in file extension: {}".format(file))
				invalid.append(os.path.basename(file))
	files_wo_ext = [base_name(x) for x in file_paths]
	file_duplicates = [x for x, y in collections.Counter(files_wo_ext).items() if y > 1]
	for f in file_duplicates:
		logger.error("{} exists multiple times. (i.e. 'foi.txt.gz' has same basename as 'foi.txt')".format(f))
		invalid.append(f)
	return invalid


def generate_randomsnps(foi_path, background, n_fois, num):
	''' Generates random SNP files in the same directory as 'foi_path' by sampling n snps randomly from the 'background'.
	'foi_paths' is used to name the randomly created snp files and where to save the newly created files
	'''
	paths = []
	out_dir = os.path.dirname(foi_path)
	if not os.path.exists(out_dir): os.mkdir(out_dir)
	for n in range(int(num)):
		rnd_snp_path = os.path.join(out_dir, "random{}_".format(n) + base_name(foi_path) + ".bed")
		# generate random snps from background
		if background.endswith('.gz'):
			command = "zcat {} | shuf -n {} | cut -f 1-3 > {}".format(background, str(n_fois), rnd_snp_path)
			out = commands.getstatusoutput(command)
		else:
			command = "shuf -n {} {}  | cut -f 1-3 > {}".format(str(n_fois), str(background), rnd_snp_path)
			out = commands.getstatusoutput(command)
		paths.append(rnd_snp_path)
	return paths


def validate_rsids(self, foi_path):
	''' Checks if the first line is contains an rsIDs. If it does, the file is sorted
	'''
	if foi_path.endswith('.gz'):
		infile = gzip.open(foi_path)
	else:
		infile = open(foi_path)
	line = infile.readline().rstrip('\n')
	infile.close()
	if line[:2] == 'rs':
		# sort the rsIDs in place
		script = "sort -k1,1 " + '-o ' + foi_path + ' ' + foi_path
		out = subprocess.Popen([script], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out.wait()
		tmp = out.stdout.read()
		tmp_er = out.stderr.read()
		if tmp_er != "":
			logger.error(tmp_er)
			raise Exception(tmp_er)
		return True
	else:
		return False



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
