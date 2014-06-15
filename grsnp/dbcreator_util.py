#!/usr/bin/env python2
import os
import subprocess
import pdb
import gzip

def remove_headers(datapath):
	''' Removes all lines that start with special characters, such as header lines.
	Replaces the datapath file with a file that does not contain header lines.
	Handles both gziped files and raw text (based on file extension)
	'''
	outputpath = datapath + ".tmp"
	if datapath.endswith('.gz'):
		infile,outfile = gzip.open(datapath),gzip.open(outputpath,'wb')
	else:
		infile,outfile = open(datapath),open(outputpath,'wb')
	while True:
		line = infile.readline().strip()
		if line == "":
			break
		if line[0].isalnum():
			outfile.write(line+"\n")
	infile.close(),outfile.close()
	os.remove(datapath)
	os.rename(outputpath,datapath)



def load_minmax(path):
	data = {}
	if not os.path.exists(path):
		return data
	score = [x for x in open(path).read().split("\n") if x != ""]
	for s in score:
		name,min_max = s.split('\t')
		data[name] = min_max
	return data


def save_minmax(data,path):
	''' Saves the dictionary of key value pairs of minmax data to a text file.
	'''
	with open(path,'wb') as writer:
		for k,v in data.items():
			writer.write("{}\t{}\n".format(k,v))


def base_name(k):
	return os.path.basename(k).split(".")[0]

def sort_convert_to_bgzip(path,outpath):
	script = "sort -k1,1 -k2,2n -k3,3n " + path +" | bgzip -c > " + outpath + ".gz.temp"
	out = subprocess.Popen([script],shell=True,stdout=subprocess.PIPE)
	out.wait()
	print 'remove ',path
	os.remove(path)# remove the .temp file extension to activate the GF		
	os.rename(outpath+".gz.temp",outpath)



def filter_by_score(gf_path_input,gf_path_output,thresh_score):
	''' Read in the gf data from gf_path_input and filter out each GF that does not
	have a score greater than the thresh_score threshold.
	gf_path_output should be WITHOUT file extension
	'''
	tmp_path = gf_path_output+'.temp'
	with open(tmp_path,"wb") as bed:
		with gzip.open(gf_path_input) as dr:
			while True:
				line = dr.readline().strip()
				if line == "":
					break
				score  = line.split('\t')[4]
				# if the score is >= to the threshold, output that GF
				if float(score) >= float(thresh_score):
					bed.write(line+"\n")
	sort_convert_to_bgzip(tmp_path,gf_path_output+'.bed.gz')

def filter_by_strand(data_dir,gf_path):
	''' Read in the gf data from gf_path and create new GFs for each strand.
	gf_path should be the full path to the gf WITHOUT file extension.
	data_dir should be the directory containing the database + grsnp_db_[score]
	'''
	sub_data_path = gf_path.replace(data_dir+"/","")
	plus_path, minus_path = os.path.join(data_dir+"_plus/",sub_data_path+'.temp'),os.path.join(data_dir+"_minus/",sub_data_path+'.temp')
	if not os.path.exists(os.path.split(plus_path)[0]): os.makedirs(os.path.split(plus_path)[0])
	if not os.path.exists(os.path.split(minus_path)[0]): os.makedirs(os.path.split(minus_path)[0])
	plus_file,minus_file = open(plus_path,'wb'),open(minus_path,'wb')
	strand_file = {"+":plus_file,'-':minus_file}
	with gzip.open(gf_path) as dr:
		while True:
			line = dr.readline().strip()
			if line == "":
				break
			strand  = line.split('\t')[5]
			# check if a valid strand exists and output to the appropriate file.
			if strand in strand_file:
				strand_file[strand].write(line+"\n")
	print plus_path, "|||" ,minus_path
	plus_file.close()
	minus_file.close()
	# remove the  strand filtered gf file if empty 
	if os.stat(plus_path).st_size==0:
		os.remove(plus_path)
	else:
		sort_convert_to_bgzip(plus_path,os.path.join(os.path.split(plus_path)[0],base_name(plus_path)+'.bed.gz'))
	if os.stat(minus_path).st_size==0:
		os.remove(minus_path)
	else:
		sort_convert_to_bgzip(minus_path,os.path.join(os.path.split(minus_path)[0],base_name(minus_path)+'.bed.gz'))


class MinMax:
	"""Used to keep track of the min and max score. 
	By Default the min and max score is None.  It is sufficient to create this class for GFs that lack a score field.
	The class will return an NA string if the score is never updated"""
	def __init__(self):
		self.max,self.min = None, None

	def update_minmax(self,n):
		
		if n == '.' or not n.isdigit():
			return
		else:
			n = float(n)
		# Assign the first value found to min and max
		if self.max == None:
			self.max = n
		if self.min == None:
			self.min = n
		# Check if we have found a new min or max
		if n < self.min:
			self.min = n
		elif n > self.max:
			self.max = n
	
	def str_minmax(self):
		n_max, n_min = 'NA','NA'
		if self.max != None:
			n_max = self.max
		if self.min != None:
			n_min = self.min
		if n_min == n_max:
			n_max, n_min = 'NA','NA'
		return "{},{}".format(n_min,n_max)

def write_line(line,path):
	with open(path, 'a') as writer:
		writer.write(line+"\n")