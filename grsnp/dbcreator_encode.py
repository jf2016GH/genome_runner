#!/usr/bin/env python2
import sys
import errno
import logging
from logging import FileHandler,StreamHandler
import os
import ftplib
import sqlite3
import string
from contextlib import closing
import subprocess
import argparse
import gzip
import re
import collections
import copy
import traceback  as trace
import pdb
from xml.sax.saxutils import quoteattr as xml_quoteattr
import BeautifulSoup
import urllib2
from grsnp.dbcreator_util import *
from time import sleep

logger = logging.getLogger('genomerunner.dbcreator')

# connection information for the ucsc ftp server
ftp_server = 'hgdownload.cse.ucsc.edu'
directory = '/goldenPath/{}/encodeDCC'
username = 'anonymous'
password = ''

download_dir = ""

class GF_ALREADY_EXISTS(Exception):
    pass

# downloads the specified file from ucsc.  Saves it with a .temp extension untill the download is complete.
def download_file(organism,gf_grp,gf_file):
	''' Downloads the gf_file from the UCSC ftp server and saves it
	in a folder with the same name as the organism.
	'''
	global download_dir
	outputpath = ''
	try:
		if os.path.exists(download_dir) == False and download_dir != '':
			logger.info( "creating directory {}".format(download_dir))
			os.makedirs(download_dir)
	except Exception, e:
		logger.warning( e)
		logger.warning("Could not create folder at {} for {}".format(download_dir,gf_file))
		return '' 
	
	try:
		outputpath = os.path.join(download_dir,gf_file)
		if not os.path.exists(outputpath):
			ftp = ftplib.FTP(ftp_server, timeout=1800) # Connection timeout 0.5h
			ftp.login(username,password)
			with open(outputpath + ".temp",'wb') as fhandle:  
				global ftp
				logger.info( 'Downloading {} from UCSC'.format(gf_file))				
				ftp.cwd(os.path.join(directory.format(organism),gf_grp))
				ftp.retrbinary('RETR ' + "{}".format(gf_file),fhandle.write)
				os.rename(outputpath+".temp",outputpath)
				logger.info( 'Finished downloading {} from UCSC'.format(gf_file))
			ftp.quit()
			# remove header lines
			remove_headers(outputpath)
		else:
			logger.info( '{} already exists, skipping download'.format(outputpath))
	except IOError as e:
		if e.errno == errno.EPIPE: # Broken pipe error handling
			logger.error('FTP download error. Restart the dbcreator. Exiting now...')
			ftp.quit()
			sys.exit(2)
	except Exception, e:
		logger.warning(e)
		logger.warning("Could not download the {} file. Names ARE case sensitive.".format(gf_file))
		return '' 
	return outputpath 


def get_gf_filepaths(organism,gf_grp):
	''' Returns the list of file in the gf_grp folder on the ucsc server.
	'''
	ftp = ftplib.FTP(ftp_server, timeout=1800) # Connection timeout 0.5h
	ftp.login(username,password)

	root_dir = os.path.join(directory.format(organism),gf_grp)
	ftp.cwd(root_dir)
	gf_paths = []
	# get all folders in the root directory
	ftp.dir(gf_paths.append)
	gf_paths = [x.split()[-1] for x in gf_paths if x[0] == '-']
	ftp.quit()
	return gf_paths



def extract_bed(organism,gf_grp,gf_file,outputpath):	
	gf_file = download_file(organism,gf_grp,gf_file)
	min_max = MinMax() # keep track of the min and max score
	if gf_file.endswith('.gz'):
		infile,outfile = gzip.open(gf_file),open(outputpath,'wb')
	else:
		infile,outfile = open(gf_file),open(outputpath,'wb')
	while True:
		line = infile.readline().rstrip('\n')
		if line == "":
			break
		line = preparebed[gf_file.replace(".gz",'').split(".")[-1]](line,min_max)+"\n"
		outfile.write(line)
	return [min_max.str_minmax(),'bed']

def _format_bed(line,min_max):
	line = line.split("\t")
	# removes special characters and the '.' used when a field is blank
	if len(line) >= 6:
		line[3] = ''.join(e for e in line[3] if e.isalnum())
		line[4] = line[4] if line[4] != "." else "0"
		line[5] = line[5] if line[5] in ["+","-"] else ""
		line = line[0:6]
		min_max.update_minmax(line[4])
	elif len(line) == 5:
		line[3] = ''.join(e for e in line[3] if e.isalnum())
		line[4] = line[4] if line[4] != "." else "0"
		line = line[0:5]
		min_max.update_minmax(line[4])
	elif len(line) == 4:
		line[3] = ''.join(e for e in line[3] if e.isalnum())
		line = line[0:4]
	return line

def _format_peak(line,min_max):
	'''Handles broadPeak and narrowPeak. Returns line of formated bed as a list of fields
	'''
	line = line.split("\t")
	# removes special characters and the '.' used when a field is blank
	line[3] = ''.join(e for e in line[3] if e.isalnum())
	line[4] = line[6] if line[6] != "." else "0" # for peak data we use the SignalValue column
	line[5] = line[5] if line[5] in ["+","-"] else ""
	line = line[0:6]
	min_max.update_minmax(line[4])
	return line

def _format_gapped_peak(line,min_max):
	'''Handles gappedPeak. Returns line of formated bed as a list of fields
	'''
	line = line.split("\t")
	# removes special characters and the '.' used when a field is blank
	line[3] = ''.join(e for e in line[3] if e.isalnum())
	line[4] = line[12] if line[12] != "." else "0" # for peak data we use the SignalValue column
	line[5] = line[5] if line[5] in ["+","-"] else ""
	line = line[0:6]
	min_max.update_minmax(line[4])
	return line

def _format_bedRNA(line,min_max):
	offset = 1
	line = line.split("\t")
	# removes special characters and the '.' used when a field is blank
	if len(line) >= 7:
		line[4] = ''.join(e for e in line[4] if e.isalnum())
		line[5] = line[5] if line[5] != "." else "0"
		line[6] = line[6] if line[6] in ["+","-"] else ""
		line = line[1:7]
		min_max.update_minmax(line[5])
	elif len(line) == 6:
		line[4] = ''.join(e for e in line[4] if e.isalnum())
		line[5] = line[5] if line[5] != "." else "0"
		line = line[1:6]
		min_max.update_minmax(line[5])
	elif len(line) == 5:
		line[4] = ''.join(e for e in line[4] if e.isalnum())
		line = line[1:5]
	return line




def get_tier(line,outputdir): # Generating paths for the ENCODE data tables using tiers, and cell types
	CELLS1 = re.compile('Gm12878|K562|H1hesc')
	CELLS2 = re.compile('A549|Cd20ro01778|Cd20ro01794|Cd20|H1neurons|Helas3|Hepg2|Huvec|Imr90|Lhcnm2|Mcf7|Monocd14ro1746|Sknsh')
	CELLS3 = re.compile('Ag04449|Ag04450|Ag09309|Ag09319|Ag10803|Aoaf|Aosmc|Be2c|Bj|Caco2|Cmk|Dnd41|Ecc1|Gm06990|Gm12801|Gm12864|Gm12865|Gm12872|Gm12873|Gm12875|Gm12891|Gm12892|Gm19239|H7es|Hac|Hae|Hah|Hasp|Hbmec|Hcfaa|Hcf|Hcm|Hcpe|Hct116|Hee|Hek293|Hffmyc|Hff|Hgf|Hipe|Hl60|Hmec|Hmf|Hmvecdblad|Hnpce|Hpae|Hpaf|Hpdlf|Hpf|Hrce|Hre|Hrpe|Hsmmfshd|Hsmmtubefshd|Hsmmt|Hsmm|Htr8|Hvmf|Jurkat|Lncap|M059j|Mcf10aes|Nb4|Nha|Nhbe|Nhdfad|Nhdfneo|Nhek|Nhlf|Nt2d1|Osteobl|Osteo|Ovcar3|Panc1|Panislets|Pfsk1|Prec|Progfib|Rpmi7951|Rptec|Saec|Skmc|Sknmc|Sknshra|T47d|Th1|Th2|U87|Werirb1|Wi38')
	m1 = CELLS1.search(line)
	m2 = CELLS2.search(line)
	m3 = CELLS3.search(line)
	if m1:
		Tier = 'Tier1'
	elif m2:
		Tier = 'Tier2'
	elif m3:
		Tier = 'Tier3'
	else:
		Tier = 'Tier3'
		with open(os.path.join(outputdir,"missing_tier.log"),'wb') as writer:
			writer.write(line)		

	return Tier


def preparebed_splitby(gf_outputdir,organism,group_name, gf_file):
	''' A function that creates separate bed files for each value in field with index 'splitby'
	gf_outputdir: the directory in which the gf_folder should be created in which the split GFs should be outputted to
	EX: /[root]/grsnp_db/[organism]/[tier]/[source]/[celltype]/

	gf_file: file name with extension i.e gfname.bed.gz
	'''
	# TODO handle the case of partially finished database
	added_features = [] 
	if not os.path.exists(gf_outputdir): os.makedirs(gf_outputdir)
	# download the GF file	
	dwnl_file = download_file(organism,	group_name,	gf_file)
	full_gf_paths = []
	# convert it to bed format
	with gzip.open(dwnl_file) as infile:
		min_max = MinMax() # keep track of the min and max score			
		file_writers = {} # {'cur_split_value': file_writer_object} A writer is created for each name field
		gf_outputdir = os.path.join(gf_outputdir,base_name(gf_file))
		if not os.path.exists(gf_outputdir):
			os.makedirs(gf_outputdir)
		while True:
			line = infile.readline().rstrip('\n')
			if line == "":
				break
			cur_gf = preparebed[gf_file.replace(".gz",'').split(".")[-1]](line,min_max)
			cur_split_value = cur_gf[3]
			# check if current TFBS already has a file writer
			outputpath = os.path.join(gf_outputdir,cur_split_value+".bed.temp")			
			if cur_split_value not in file_writers:
				o_dir = os.path.dirname(outputpath)
				new_path = os.path.join(o_dir,''.join(e for e in base_name(outputpath) if e.isalnum() or e=='.' or e=='_')) + ".bed.gz"
				if os.path.exists(new_path):
					raise GF_ALREADY_EXISTS("{} already exists. Not going to overwrite it.".format(new_path))
				file_writers[cur_split_value] = open(outputpath,'wb')
				full_gf_paths.append(outputpath)
			file_writers[cur_split_value].write("\t".join(cur_gf)+"\n")
		#close all open files
		for k in file_writers.keys():
			file_writers[k].close()
		# convert all created files into gzip
		converted_paths = []
		for f_path in full_gf_paths:
			sort_convert_to_bgzip(f_path,new_path)
			converted_paths.append(new_path)
		return [min_max.str_minmax(),converted_paths]

def preparebed(gf_outputdir, organism, group_name, gf_file):
	''' Converts the file to the correct bed format, sorts it, and gzips it. Returns the min_max stats, 
	outputdir: the directory in which the gf_folder should be created in which the split GFs should be outputted to
	EX: /[root]/grsnp_db/[organism]/[tier]/[source]/[celltype]/

	gf_file: file name with extension i.e gfname.bed.gz
	'''
	added_features = [] 
	if not os.path.exists(outputdir): os.makedirs(outputdir)
	# download the GF file	
	dwnl_file = download_file(organism,	group_name,	gf_file)
	f_path = os.path.join(outputdir,base_name(gf_file)+".bed.temp")
	o_dir = os.path.dirname(f_path)
	new_path = os.path.join(o_dir,''.join(e for e in base_name(f_path) if e.isalnum() or e=='.' or e=='_')) + ".bed.gz"
	if os.path.exists(new_path) == True:
		raise GF_ALREADY_EXISTS("{} already exists. Not going to overwrite it.".format(new_path))
	# convert it to bed format
	with gzip.open(dwnl_file) as infile:
		min_max = MinMax() # keep track of the min and max score
		with open(f_path,'wb') as writer:
			while True:
				line = infile.readline().rstrip('\n')
				if line == "":
					break
				cur_gf = preparebed[gf_file.replace(".gz",'').split(".")[-1]](line,min_max)
				writer.write("\t".join(cur_gf)+"\n")			
			# convert all created files into gzip
			sort_convert_to_bgzip(f_path,new_path)
	return [min_max.str_minmax(),[new_path]]





def _get_celltype_source(f_name, padding):
	''' Extracts the source and cell type from the genomic feature file name
	'''
	f_name.lstrip(padding)
	categories = re.findall('[A-Z][^A-Z]*', f_name)
	return {"source": categories[0], "cell": categories[1]}

def _get_gf_directory(outputdir,gf_group,gf_name):
	''' Returns the output_dir of the gf.
	'''
	# Dictates how much of the filename to strip off before searching for cell type/source etc.
	padding = {
	 "wgEncodeAwgTfbsUniform": "wgEncodeAwgTfbs",
	 "wgEncodeBroadHmm": "wgEncodeBroad",
	 "wgEncodeBroadHistone": "wgEncodeBroadHistone",
	 "wgEncodeUwHistone": "wgEncodeUwHistone",
	 "wgEncodeSydhHistone": "wgEncodeSydhHistone",
	 "wgEncodeAwgDnaseUniform": "wgEncodeAwgDnaseUniform"
	}

	# Dictates the structure of the directory
	dirstruture = {
	 "wgEncodeAwgTfbsUniform": ['tier','cell','source'],
	 'wgEncodeBroadHmm': ['tier','cell'],
	 "wgEncodeRegTfbsClustered": [],
	 "wgEncodeBroadHistone": ['tier','cell','source'],
	 "wgEncodeUwHistone": ['tier','cell','source'],
	 "wgEncodeSydhHistone": ['tier','cell','source'],
	 "wgEncodeAwgDnaseUniform": ['tier','cell','source']
	}

	root_folder = {
	"wgEncodeAwgTfbsUniform": "ENCODE",
	 "wgEncodeRegTfbsClustered": "ENCODE",
	 "wgEncodeBroadHmm": "ENCODE",
	 "wgEncodeBroadHistone": "ENCODE",
	 "wgEncodeUwHistone": "ENCODE",
	 "wgEncodeSydhHistone": "ENCODE",
	 "wgEncodeAwgDnaseUniform": "ENCODE"
	}

	dir_structure = dirstruture[gf_group]
	gf_directory = [root_folder[gf_group]]
	for folder in dirstruture[gf_group]:
		if folder == 'tier':
			gf_directory.append(get_tier(gf_name,outputdir))
		elif folder == 'cell':
			gf_directory.append(_get_celltype_source(gf_name,padding[gf_group])['cell'])
		elif folder == 'source':
			gf_directory.append(_get_celltype_source(gf_name,padding[gf_group])['source'])
	gf_directory = "/".join(gf_directory)
	return os.path.join(outputdir,gf_directory)

def create_feature_set(data_dir,organism,gf_group,pct_score=None,max_install = None):
	outputdir = os.path.join(data_dir,organism)
	min_max_path = os.path.join(outputdir,'minmax.txt')
	added_features = [] 
	outpath = ""
	prog, num = 0,len(gfs)
	summary_path = os.path.join(outputdir,"summary.log")
	if not os.path.exists(outputdir): os.makedirs(outputdir)
	open(summary_path,'wb')


	min_max_scores = load_minmax(min_max_path)
	grp_count,gf_file_paths = 0,[]
	if "gf_files" in prepare_type[gf_group].keys():
		gf_file_paths = prepare_type[gf_group]["gf_files"]
	else:
		gf_file_paths = get_gf_filepaths(organism,gf_group)
	for gf_file in gf_file_paths:
		# check if gf_type is supported
		if gf_file.replace(".gz",'').split(".")[-1] not in preparebed.keys():
			continue
		# limit the number of GFs to install per group	
		if max_install != None and grp_count >= max_install:
			break
		try:
			gf_type = ""
			gf_outputdir = _get_gf_directory(outputdir,gf_group,base_name(gf_file))
			if not os.path.exists(gf_outputdir):
				os.makedirs(gf_outputdir)
			# removes the .temp file, to prevent duplicate data from being written
			if os.path.exists(outpath+".temp"):
				os.remove(outpath+".temp")
			# converts the ucsc data into proper bed format
			logger.info( "Converting into proper bed format: {}".format(gf_file))
			try:
				[minmax_score, gf_paths] = prepare_type[gf_group]["prep_method"](gf_outputdir,organism,gf_group,gf_file)
			except GF_ALREADY_EXISTS:
				logger.info( "{} already exists as or .gz, skipping extraction".format(outpath.replace(".gz","")))
				continue
			# output minmax stats
			for f in gf_paths:
				min_max_scores[f] = minmax_score
			save_minmax(min_max_scores,min_max_path)

			# cannot detect type, skip
			if gf_type == "failed":
				write_line("\t".join([gf_file,"Not supported","None"]),summary_path)
				logger.warning( "Unable to convert {} into bed".format(gf_file))
				continue
			added_features.append(outpath)
			grp_count += 1			
			write_line("\t".join([gf_file,gf_type,"None"]),summary_path)

		except Exception, e:
			write_line("\t".join([gf_file,"Failed",str(e)]),summary_path)
			exc = trace.format_exc()
			logger.warning( "Unable to convert {} into bed".format(gf_file))
			logger.warning(exc)
			prog += 1
			continue	

	prog += 1
	# cleanup the temporary files
	if os.path.exists(outpath + ".temp"): os.remove(outpath+".temp")

	# logger.info( "The following types are not supported (includes all 'big' file types):\n " + str(notsuptypes))
	# logger.info("The following features were added to the database: \n{}".format(added_features))
	return "Created UCSC database"

def update_progress(progress):
    print '\r[{0}] {1}%'.format('#'*(progress/10), progress),


# each genomic feature group must have an entry in this dictionary
# f_names is optional and can be left out.  It limits the download of the gf group to only the files listed
prepare_type = {
	"wgEncodeRegTfbsClustered": {"prep_method":  preparebed_splitby, "gf_files": ["wgEncodeRegTfbsClusteredV3.bed.gz"]},
	"wgEncodeBroadHmm": {"prep_method":  preparebed_splitby},
	"wgEncodeBroadHistone": {"prep_method":  preparebed},
	"wgEncodeUwHistone": {"prep_method":  preparebed},
	"wgEncodeSydhHistone": {"prep_method":  preparebed},
	"wgEncodeAwgDnaseUniform": {"prep_method":  preparebed}
}


preparebed = {
	"bed":_format_bed,
	"bed9":_format_bed,
	"bedRrbs":_format_bed,
	"peptideMapping":_format_bed,
	"broadPeak":_format_peak,
	"narrowPeak":_format_peak,
	"gappedPeak": _format_gapped_peak,
	"bedRnaElements":_format_bed

}

if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog="python -m grsnp.dbcreator", description='Creates the GenomeRunner SNP Database. Example: python -m grsnp.dbcreator -d /home/username/grsnp_db/ -g mm9', epilog='IMPORTANT: Execute DBCreator from the database folder, e.g., /home/username/grsnp_db/. Downloaded files from UCSC are placed in ./downloads database created in ./grsnp_db.')
	parser.add_argument("--data_dir" , "-d", nargs="?", help="Set the directory where the database to be created. Use absolute path. Example: /home/username/grsnp_db/. Required", required=True)
	parser.add_argument('--organism','-g', nargs="?", help="The UCSC code of the organism to use for the database creation. Default: hg19 (human). Required", default="hg19")
	parser.add_argument('--featuregroups','-f', nargs="?", help='The names of the specific genomic feature groups to download.  List available for hg19 at ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/', default="")
	parser.add_argument('--galaxy', help="Create the xml files needed for Galaxy. Outputted to the current working directory.", action="store_true")
	parser.add_argument('--score', '-s', help="Commas separated list of score percentiles.", nargs='?',default="")
	#parser.add_argument('--filteronly','-o', help="Only filter by score and strand. Skips downloading and installing new GFs.", action="store_true")
	parser.add_argument('--max','-m', nargs="?", help="Limit the number of features to be created within each group.",type=int, default=None)



	args = vars(parser.parse_args())

	hdlr = logging.FileHandler(os.path.join(args['data_dir'], 'genomerunner_dbcreator.log'))
	hdlr_std = StreamHandler()
	formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
	hdlr.setFormatter(formatter)
	logger.addHandler(hdlr)
	logger.addHandler(hdlr_std)
	logger.setLevel(logging.INFO)

	if not args["data_dir"]:
		print "ERROR: --data_dir is required"
		sys.exit()

	if args['score'] == "":
		args['score'] = "25,50,75"
	args['score'] = set(args['score'].split(',')) # remove duplicate scores
	data_dir=os.path.join(args["data_dir"],'grsnp_db')


	if args['organism'] is not None: # Only organism is specified. Download all organism-specific features
		global download_dir, prepare_type
		download_dir = os.path.join(args["data_dir"],"downloads",args['organism'])
		gfs = args["featuregroups"].split(",")
		for grp in prepare_type.keys():
			create_feature_set(data_dir,args['organism'],grp,None,5)
	else:
		print "ERROR: Requires UCSC organism code.  Use --help for more information"
		sys.exit()

	# load score from minmax.txt file created earlier
	minmax = load_minmax(os.path.join(data_dir,args['organism'],"minmax.txt"))		


	### Second Step: Create subdirectories for score and filter data by score percentile
	# create sub directories for score percentiles and populate with score-filtered GF data
	# gather all directories (groups) in the database
	sys.exit(0)
	print "Filtering GFs by strand and score..."
	orgdir = os.path.join(data_dir,args['organism'])
	dirs = [name for name in os.listdir(orgdir)
		if os.path.isdir(os.path.join(orgdir, name))]
	for d in dirs:
		# gather all paths
		gfs = []
		for base, tmp, files in os.walk(os.path.join(orgdir,d)):
				gfs += [os.path.join(base,f) for f 	 
					in files if f.endswith(('.gz'))]
		for gf_path in gfs: 
			print "Filtering {} ...".format(gf_path)
			# filter the original GF by strand
			filter_by_strand(data_dir,gf_path)
			for pct_score in args['score']:		
				[score_min,score_max] = minmax[base_name(gf_path)].split(",")
				# calculate threshold score
				if score_min == 'NA':
					continue
				score_min,score_max = float(score_min),float(score_max) 
				thresh_score = score_min + (score_max-score_min)*float(pct_score)/100
				logger.info("MinMax stats for {}: Min={}, Max={}, {} pct_thresh={}".format(base_name(gf_path), score_min,score_max,pct_score,thresh_score))
				# is this safe? It searches /dirpath/grsnp_db/subdirs/gf.txt and replaces /grsnp_db/ with /grsnp_db_[score]/
				gf_scorepath_out =gf_path.replace('/grsnp_db/','/grsnp_db_{}/'.format(pct_score))
				if not os.path.exists(os.path.split(gf_scorepath_out)[0]):
					os.makedirs(os.path.split(gf_scorepath_out)[0])
				gf_path_out_woext = os.path.join(os.path.split(gf_scorepath_out)[0],base_name(gf_scorepath_out))
				# filter by score
				filter_by_score(gf_path, gf_path_out_woext,thresh_score)
				# filter the score filtered GF by strand				
				filter_by_strand(data_dir+"_{}".format(pct_score),gf_scorepath_out)


	root_dir = os.path.dirname(os.path.realpath(__file__))
	readme = open(os.path.join(root_dir,"grsnp_db_readme.txt")).read()
	with open("grsnp_db_readme.txt","wb") as writer:
		writer.write(readme)
	print "FINISHED: Downloaded files from UCSC are placed in {}.  Database created in {}".format(os.path.join(args["data_dir"],"downloads"),os.path.join(args["data_dir"],"grsnp_db`"))




