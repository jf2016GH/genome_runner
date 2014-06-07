#!/usr/bin/env python2
import os
import subprocess


def remove_headers(datapath):
	''' Removes all lines that start with special characters, such as header lines.
	Replaces the datapath file with a file that does not contain header lines.
	Handles both gziped files and raw text (based on file extension)
	'''
	outputpath = datapath + ".tmp"
	if datapath[:-3] == '.gz':
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
	os.remove(path)# remove the .temp file extension to activate the GF		
	os.rename(outpath+".gz.temp",outpath)

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