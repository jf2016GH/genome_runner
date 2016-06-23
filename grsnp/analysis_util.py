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
import datetime
DEBUG = False


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


def _zip_run_files(outdir, job_id=""):
	'''
	File paths of FOIs and GFs as a list. Gathers all the files together in one zipped file
	'''
	# zip annotation result folder if it exists
	anno_dir = os.path.join(outdir, 'annotations')
	if os.path.exists(anno_dir):
		anno_zip_path = os.path.join(outdir, 'annotations.zip')
		if os.path.exists(anno_zip_path):
			os.remove(anno_zip_path)
		z_ano = zipfile.ZipFile(anno_zip_path, 'a')
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
	tar_path = os.path.join(outdir, 'GR_{}.tar'.format(job_id))
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


def validate_rsids(foi_path):
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


def _write_head(content, outpath):
	f = front_appender(outpath)
	f.write(content)
	f.close()

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





def _chunks(l, n):
	n = max(1, n)
	return [l[i:i + n] for i in range(0, len(l), n)]


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
