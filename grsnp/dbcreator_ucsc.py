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

try:
	import xml.etree.cElementTree as ET
except ImportError:
	import xml.etree.ElementTree as ET
# connection information for the ucsc ftp server
ftp_server = 'hgdownload.cse.ucsc.edu'
directory = '/goldenPath/{}/database'
username = 'anonymous'
password = ''

logger = logging.getLogger('genomerunner.dbcreator')

ftp = ""
illegal_chars = ['=',':']


			
# downloads the specified file from ucsc.  Saves it with a .temp extension untill the download is complete.
def download_ucsc_file(organism,filename,downloaddir,remove_headers=True):
	''' Downloads the filename from the UCSC ftp server and saves it
	in a folder with the same name as the organism.
	'''
	outputpath = ''
	if downloaddir != None and downloaddir != '':
		outputdir = os.path.join(downloaddir,organism)
	else:
		outputdir = organism
	try:
		if os.path.exists(outputdir) == False and outputdir != '':
			logger.info( "creating directory {}".format(outputdir))
			os.makedirs(outputdir)
	except Exception, e:
		logger.warning( e)
		logger.warning("Could not create folder at {} for {}".format(outputdir,filename))
		return '' 
	
	try:
		outputpath = os.path.join(outputdir,filename)
		if not os.path.exists(outputpath):
			with open(outputpath + ".temp",'wb') as fhandle:  
				global ftp
				logger.info( 'Downloading {} from UCSC'.format(filename))				
				ftp.cwd(directory.format(organism))
				ftp.retrbinary('RETR ' + "{}".format(filename),fhandle.write)
				os.rename(outputpath+".temp",outputpath)
				logger.info( 'Finished downloading {} from UCSC'.format(filename))			
			# remove header lines
			if remove_headers:
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
		logger.warning("Could not download the {} sql file. Names ARE case sensitive.".format(filename))
		return '' 

	return outputpath 


def download_trackdb(organism,outputdir):
	''' Downloads the trackdb.sql and trackDb.txt.gz from the UCSC ftp server and saves it in a folder with the same name as the organism.
		Returns the path of the downloaded .sql file
	'''
	sqloutputpath = download_ucsc_file(organism,"trackDb.sql",outputdir,False)
	dataoutpath = download_ucsc_file(organism,"trackDb.txt.gz",outputdir,False)
	' replace all of the \\\n characters in the html column with <br />'
	text = gzip.open(dataoutpath).read()
	with gzip.open(dataoutpath,'wb') as sw:
		sw.write(text.replace('\\\n','<br />').replace('\\\t','     '))
	return sqloutputpath


	

def _check_cols(colnames,colstoextract):
	'''checks if all the columns to extract
	actually exist in the ucsc table
	'''
	for c in colstoextract:
		if not c in colnames:
			return False
	return True


def extract_bed6(outputpath,datapath,colnames):
	colstoextract,mm = ['chrom','chromStart','chromEnd','name','score','strand'],MinMax()
	# Checks if all of the columns exist in the table.  If not extract_bed5 is tried i
	if _check_cols(colnames,colstoextract):
		logger.info( "Outpath is: {}".format(outputpath))
		with gzip.open(datapath) as dr:
			with open(outputpath,"wb") as bed:
				while True:
					line = dr.readline().strip('\r').rstrip('\n')
					if line == "":
						break
					r  = dict(zip(colnames,line.split('\t')))
					row = []
					row = [r["chrom"],r["chromStart"],r["chromEnd"],''.join(e for e in r["name"] if e.isalnum()),r["score"] if r["score"] != "." else "0",r["strand"] if r["strand"] in ["+","-"] else ""]# Can't use strand as "."
					bed.write("\t".join(map(str,row))+"\n")
					# check if new min or max score found
					mm.update_minmax(r['score'])
		return [mm.str_minmax(),"bed 6"]
	else:
		logger.warning("Nonstandard bed6, attempting extraction as bed5")
		return extract_bed5(outputpath,datapath,colnames)

def extract_bed5(outputpath,datapath,colnames):
	colstoextract,mm = ['chrom','chromStart','chromEnd','name','score'],MinMax()
	if _check_cols(colnames,colstoextract):
		with open(outputpath,"wb") as bed:
			with gzip.open(datapath) as dr:
				while True:
					line = dr.readline().rstrip('\r').rstrip('\n')
					if line == "":
						break
					r  = dict(zip(colnames,line.split('\t')))
					row = [r["chrom"],r["chromStart"],r["chromEnd"],''.join(e for e in r["name"] if e.isalnum()),r["score"] if r["score"] != "." else "0"] # Can't use strand as "."
					bed.write("\t".join(map(str,row))+"\n")
					# check if new min or max score found
					mm.update_minmax(r['score'])
		return [mm.str_minmax(),"bed 5"]
	else:
		logger.warning("Nonstandard bed5, attempting extraction as bed4")
		return extract_bed4(outputpath,datapath,colnames)
	
def extract_bed4(outputpath,datapath,colnames):
	colstoextract,mm = ['chrom','chromStart','chromEnd','name'],MinMax()
	if _check_cols(colnames,colstoextract):
		with open(outputpath,"wb") as bed:
			with gzip.open(datapath) as dr:
				while True:
					line = dr.readline().rstrip('\r').rstrip('\n')
					if line == "":
						break
					r  = dict(zip(colnames,line.split('\t')))
					row = [r["chrom"],r["chromStart"],r["chromEnd"],''.join(e for e in r["name"] if e.isalnum()).replace(": ",""),"0"] # Can't use strand as ".". Replace ": " is needed for cpgIslandExt
					bed.write("\t".join(map(str,row))+"\n")
		return [mm.str_minmax(),"bed 4"]
	else:
		logger.warning("Nonstandard bed4, attempting extraction as bed3")
		return extract_bed3(outputpath,datapath,colnames)

def extract_bed3(outputpath,datapath,colnames):
	colstoextract,mm = ['chrom','chromStart','chromEnd'],MinMax()
	if _check_cols(colnames, colstoextract):
		with open(outputpath,"wb") as bed:
			with gzip.open(datapath) as dr:
				while True:
					line = dr.readline().strip('\r').rstrip('\n')
					if line == "":
						break
					r  = dict(zip(colnames,line.split('\t')))
					row = [r["chrom"],r["chromStart"],r["chromEnd"],"","0"] # Can't use strand as "."
					bed.write("\t".join(map(str,row))+"\n")
		return [mm.str_minmax(),"bed 3"]
	else:
		logger.warning("bed3 failed.")
		return [mm.str_minmax(),"failed"]

def extract_psl(outputpath,datapath,colnames):
	colstoextract,mm = ['tName','tStart','tEnd','qName','qSize','strand'],MinMax()
	# Checks if all of the columns exist in the table.  If not
	if _check_cols(colnames,colstoextract):
		logger.info( "Outpath is: {}".format(outputpath))
		with gzip.open(datapath) as dr:
			with open(outputpath,"wb") as bed:
				while True:
						line = dr.readline().strip('\r').rstrip('\n')
						if line == "":
							break
						r  = dict(zip(colnames,line.split('\t')))
						row = []
						row = [r["tName"],r["tStart"],r["tEnd"],''.join(e for e in r["qName"] if e.isalnum()),r["qSize"] if r["qSize"] != "." else "0",r["strand"] if r["strand"] in ["+","-"] else ""]# Can't use strand as "."
						bed.write("\t".join(map(str,row))+"\n")
		return [mm.str_minmax(),"psl"]
	else:
		return [mm.str_minmax(),"failed"]
		
def extract_genepred(outputpath,datapath,colnames):
	colstoextract,mm = ['chrom','txStart','txEnd','name','strand'],MinMax()
	exonpath = outputpath.split(".")[0]+"_exon"
	if _check_cols(colnames,colstoextract):
		# removes the .temp file of the exon, to prevent duplicate data from being written
		if os.path.exists(exonpath+".temp"): 
			os.remove(exonpath+".temp")
		with gzip.open(datapath) as dr:
			with open(outputpath,"wb") as bed:
				with open(exonpath+".temp","wb") as exonbed:				
					while True:
						line = dr.readline().rstrip('\r').rstrip('\n')
						if line == "":
							break
						r = dict(zip(colnames,line.split('\t')))
						# extract the gene data inserts a blank for score
						row = [r['chrom'],r['txStart'],r['txEnd'],''.join(e for e in r['name'] if e.isalnum()),'0',r['strand']]
						bed.write("\t".join(map(str,row))+"\n")
						
						# extract the exon data
						for (s,e) in zip(r["exonStarts"].split(","),r["exonEnds"].split(",")):
							if s != '':
								rowexon = [r['chrom'],s,e,''.join(e for e in r['name'] if e.isalnum()),'0',r['strand']]
								exonbed.write("\t".join(map(str,rowexon))+"\n")
		# sort the exon file and convert to bgzip format
		sort_convert_to_bgzip(exonpath+".temp",exonpath + ".bed.gz")
		return [mm.str_minmax(),"genepred"]
	else:
		return [mm.str_minmax(),"failed"]

def extract_peak(outputpath,datapath,colnames):
	colstoextract,mm = ['chrom','chromStart','chromEnd','name','signalValue','strand'],MinMax()
	# Checks if all of the columns exist in the table.  If not extract_bed5 is tried i
	if _check_cols(colnames,colstoextract):
		logger.info( "Outpath is: {}".format(outputpath))
		with gzip.open(datapath) as dr:
			with open(outputpath,"wb") as bed:
				while True:
					line = dr.readline().strip('\r').rstrip('\n')
					if line == "":
						break
					r  = dict(zip(colnames,line.split('\t')))
					row = []
					row = [r["chrom"],r["chromStart"],r["chromEnd"],''.join(e for e in r["name"] if e.isalnum()),r["signalValue"] if r["signalValue"] != "." else "0",r["strand"] if r["strand"] in ["+","-"] else ""]# Can't use strand as "."
					bed.write("\t".join(map(str,row))+"\n")
					# check if new min or max score found
					mm.update_minmax(r['signalValue'])
		return [mm.str_minmax(),"peak"]
	else:
		return [mm.str_minmax(),"failed"]


def guess_extract_type(outpath,datapath,colnames):
	[minmax,gf_type] = extract_peak(outpath,datapath,colnames)
	if gf_type == "failed": [minmax,gf_type] = extract_bed6(outpath,datapath,colnames)
	if gf_type == "failed": [minmax,gf_type] = extract_psl(outpath,datapath,colnames)
	if gf_type == "failed": [minmax,gf_type] = extract_genepred(outpath,datapath,colnames)
	return [minmax,gf_type]



def extract_rmsk(outputpath,datapath,colnames):
	colstoextract,mm = ['genoName','genoStart','genoEnd','repClass', 'strand','swScore'],MinMax()
	with open(outputpath,"wb") as bed:
		with gzip.open(datapath) as dr:
			while True:
				line = dr.readline().strip('\r').rstrip('\n')
				if line == "":
					break
				r = dict(zip(colnames,line.split('\t')))
				row = [r["genoName"],r["genoStart"],r["genoEnd"],''.join(e for e in r["repClass"] if e.isalnum()),r["swScore"],r["strand"]]
				bed.write("\t".join(map(str,row))+"\n")
				mm.update_minmax(r['swScore'])
	return [mm.str_minmax(),"rmsk"]

def get_column_names(sqlfilepath):
	''' extracts the column names from the .sql file and returns them
	'''
	tdbsql = open(sqlfilepath).read()
	# creates a tuple containing the column names
	tbdcolumns = re.findall("\n\s\s`*(.+?)`*\s", tdbsql, re.DOTALL)
	tbdcolumns = [c for c in tbdcolumns if not c=="KEY"]
	return tbdcolumns	


# The different file types that can be extracted.
# To add a new type add a new entry into this dictionary along with 
# the name of the function that should be used to extract the data.
preparebed = {"bed 6" : extract_bed6,
				"broadPeak": extract_peak,
				"narrowPeak": extract_peak,
				"bed 6 +" : extract_bed6,
				"bed 12 +": extract_bed6,
				"bed 12 .": extract_bed6,
				"bed 12": extract_bed6,
				"bed 10": extract_bed6,
				"bed 9 +": extract_bed6,
				"bed 9 .": extract_bed6,
				"bed 9": extract_bed6,
				"bed 8 +": extract_bed6,
				"bed 8 .": extract_bed6,
				"bed 6 .": extract_bed6,
				"bed 5 +": extract_bed5,
				"bed 5 .": extract_bed5,
				"bed 5": extract_bed5,
				"bed5FloatScore": extract_bed5,
				"bed 4 +": extract_bed4,
				"bed 4 .": extract_bed4,
				"bed 4": extract_bed4,
				"bed 3 +": extract_bed3,
				"bed 3 .": extract_bed3,
				"bed 3" : extract_bed3,
				"genePred xenoRefPep xenoRefMrna": extract_genepred,
				"genePred vegaPep": extract_genepred,
				"genePred sgpPep": extract_genepred,
				"genePred refPep refMrna": extract_genepred,
				"genePred nscanPep": extract_genepred,
				"genePred knownGenePep knownGeneMrna": extract_genepred,
				"genePred genscanPep": extract_genepred,
				"genePred geneidPep": extract_genepred,
				"genePred ensPep": extract_genepred,
				"genePred acemblyPep acemblyMrn": extract_genepred,
				"genePred acemblyPep acemblyMrna": extract_genepred,
				"genePred" : extract_genepred,
				"psl" : extract_psl,
				"psl ." : extract_psl,
				"psl est" : extract_psl,
				"psl protein" : extract_psl,
				"psl xeno" : extract_psl,
				"rmsk" : extract_rmsk,
				"factorSource" : extract_bed6,
				None: guess_extract_type}
				
numdownloaded = collections.defaultdict(int)

def encodePath(line): # Generating paths for the ENCODE data tables using groups, tiers, and cell types
	ENCODE = re.compile('AffyRnaChipFiltTransfrags|BroadHistone|BroadHmm|GisChiaPet|GisRnaPet|HaibMethyl450|HaibGenotype|HaibMethylRrbs|HaibTfbs|OpenChromSynth|RikenCage|SunyAlbanyGeneSt|SunyAlbanyTiling|SunyRipSeq|SunySwitchgear|UmassDekker5C|UwAffyExonArray|UwDgf|UwDnase|UwHistone|UwRepliSeq|UwTfbs|CshlLongRnaSeq|CshlShortRnaSeq|LicrHistone|LicrTfbs|PsuHistone|PsuTfbs')
	CELLS1 = re.compile('Gm12878|K562|H1hesc')
	CELLS2 = re.compile('A549|Cd20ro01778|Cd20ro01794|Cd20|H1neurons|Helas3|Hepg2|Huvec|Imr90|Lhcnm2|Mcf7|Monocd14ro1746|Sknsh')
	CELLS3 = re.compile('Ag04449|Ag04450|Ag09309|Ag09319|Ag10803|Aoaf|Aosmc|Be2c|Bj|Caco2|Cmk|Dnd41|Ecc1|Gm06990|Gm12801|Gm12864|Gm12865|Gm12872|Gm12873|Gm12875|Gm12891|Gm12892|Gm19239|H7es|Hac|Hae|Hah|Hasp|Hbmec|Hcfaa|Hcf|Hcm|Hcpe|Hct116|Hee|Hek293|Hffmyc|Hff|Hgf|Hipe|Hl60|Hmec|Hmf|Hmvecdblad|Hnpce|Hpae|Hpaf|Hpdlf|Hpf|Hrce|Hre|Hrpe|Hsmmfshd|Hsmmtubefshd|Hsmmt|Hsmm|Htr8|Hvmf|Jurkat|Lncap|M059j|Mcf10aes|Nb4|Nha|Nhbe|Nhdfad|Nhdfneo|Nhek|Nhlf|Nt2d1|Osteobl|Osteo|Ovcar3|Panc1|Panislets|Pfsk1|Prec|Progfib|Rpmi7951|Rptec|Saec|Skmc|Sknmc|Sknshra|T47d|Th1|Th2|U87|Werirb1|Wi38')
	n = ENCODE.search(line) 
	m1 = CELLS1.search(line)
	m2 = CELLS2.search(line)
	m3 = CELLS3.search(line)
	if n:
		grp = n.group()
	else:
		grp = 'Special'
	if m1:
		Tier = 'Tier1'
		Cell = m1.group()
	elif m2:
		Tier = 'Tier2'
		Cell = m2.group()
	elif m3:
		Tier = 'Tier3'
		Cell = m3.group()
	else:
		Tier = 'Tier3'
		Cell = ''
	return os.path.join('ENCODE', grp, Tier, Cell, line.strip())		

def create_feature_set(trackdbpath,organism,max_install,gfs=[],pct_score=None):
	outputdir = os.path.dirname(trackdbpath)
	download_dir = os.path.join(os.path.split(os.path.split(outputdir)[0])[0],"downloads")
	min_max_path = os.path.join(outputdir,'minmax.txt')
	trackdb = load_tabledata_dumpfiles(os.path.splitext(trackdbpath)[0])
	added_features = [] 
	notsuptypes, outpath = set([]),""
	prog, num = 0,len(gfs)
	summary_path = os.path.join(outputdir,"summary.log")
	gfs = [x for x in gfs if x != ""]
	# get list of GFs files on server
	if len(gfs) == 0:
		html = urllib2.urlopen('http://hgdownload.cse.ucsc.edu/goldenPath/{}/database/'.format(organism)).read()
		soup = BeautifulSoup.BeautifulSoup(html)
		gfs = [x.get('href')[:-4] for x in soup.findAll('a') if 'sql' in x.get('href')]

	min_max_scores,num = load_minmax(min_max_path),len(gfs)
	for gf_name in gfs:
		# check if GF is listed in trackdb
		row = [x for x in trackdb if x['tableName'] == gf_name]
		# if in trackdb, process using type provided
		if len(row) != 0:
			#continue
			row = row[0] # convert from list to dictionary
		else:
			#if 'chain' in gf_name: continue
			row = {'type': None,'grp':'unsorted','tableName': gf_name,'longLabel':'None'}
			logger.info( 'Processing files {} of {}'.format(prog,num))
		# if table name has Raw, in it, it is likely raw reads, skip
		if _exclude(gf_name):
			write_line("\t".join([gf_name,"Raw signal",row["longLabel"]]),summary_path)
			logger.info('Raw signal detected: {}, skipping...'.format(gf_name))
			prog +=1
			continue
		if row['type'] in preparebed:
			# this line limits the number of GFs to download, note that this check only occurs for GFs in tracdb
			if numdownloaded[str(row["type"])] <= max_install or max_install == None:
				sqlpath = download_ucsc_file(organism,row["tableName"] + ".sql",download_dir,False)
				download_ucsc_file(organism,row["tableName"] + ".txt.gz",download_dir,True)
				if sqlpath != '':
					try:
						gf_type = ""
						if row["tableName"].startswith("wgEncode"):
							outpath = os.path.join(outputdir, encodePath(row["tableName"]))
						else:
							outpath = os.path.join(outputdir,row["grp"],row["tableName"]) # ,'Tier' + row["visibility"]
						if not os.path.exists(os.path.dirname(outpath)):
							os.makedirs(os.path.dirname(outpath))
						if os.path.exists(outpath + ".bed.gz") == False:
							# removes the .temp file, to prevent duplicate data from being written
							if os.path.exists(outpath+".temp"):
								os.remove(outpath+".temp")
							# converts the ucsc data into propery bed format
							logger.info( "Converting into proper bed format: {}".format(os.path.splitext(sqlpath)[0]))
							[minmax_score,gf_type] = preparebed[row["type"]](outpath+".temp",os.path.splitext(sqlpath)[0]+".txt.gz",get_column_names(os.path.splitext(sqlpath)[0]+".sql"))
							# cannot detect type, skip
							if gf_type == "failed":
								write_line("\t".join([gf_name,"Not supported","None"]),summary_path)
								logger.warning( "Unable to convert {} into bed".format(gf_name))
								prog += 1
								continue
							# output minmax stats
							min_max_scores[row['tableName']] = minmax_score
							save_minmax(min_max_scores,min_max_path)
							# sort the file and convert to bgzip format
							o_dir = os.path.dirname(outpath)
							new_path = os.path.join(o_dir,''.join(e for e in os.path.basename(outpath) if e.isalnum() or e=='.' or e=='_')) + ".bed.gz"
							sort_convert_to_bgzip(outpath+".temp",new_path)
							added_features.append(outpath)
						else:
							logger.info( "{} already exists as or .gz, skipping extraction".format(outpath.replace(".gz","")))
						write_line("\t".join([gf_name,gf_type,row["longLabel"]]),summary_path)
						numdownloaded[str(row["type"])] += 1

					except Exception, e:
						write_line("\t".join([gf_name,"Failed",str(e)]),summary_path)
						exc = trace.format_exc()
						logger.warning( "Unable to convert {} into bed".format(row["tableName"]))
						logger.warning(exc)
						prog += 1
						continue
		else:
			write_line("\t".join([gf_name,"{} Not supported".format(row["type"]),row["longLabel"]]),summary_path)
			if 'big' not in row['type']:
				notsuptypes.add(row['type'])	

	prog += 1
	# cleanup the temporary files
	if os.path.exists(outpath + ".temp"): os.remove(outpath+".temp")

	# logger.info( "The following types are not supported (includes all 'big' file types):\n " + str(notsuptypes))
	# logger.info("The following features were added to the database: \n{}".format(added_features))
	logger.info("A count of features added by type: ")
	for k,d in numdownloaded.iteritems():
		logger.info( k + ":" + str(d))
	return "created database"
	

def write_line(line,path):
	with open(path, 'a') as writer:
		writer.write(line+"\n")


def _gettype(feature,trackdb):
	'''Returns the type of feature from the trackdb'''
	for x in trackdb:
		if x['tableName'] == feature:
			return x['type']

def _get_info(feature,trackdb):
	''' Returns the row in trackdb that contains the feature information'''
	for t in trackdb:
		if feature == t['tableName']:
			return t
	return False



def load_tabledata_dumpfiles(datapath):
	''' Loads the table data into memory from the sql file and the .txt.gz file
	data path should be WITHOUT extension (example. home/trackDb
	'''
	colnames = get_column_names(datapath+".sql")
	data = list()	
	with gzip.open(datapath+'.txt.gz') as fhandle:
		while True:
			line = fhandle.readline()
			if line == "":
				break
			row  = dict(zip(colnames,line.split('\t')))
			data.append(row)
	return data


def create_galaxy_xml_files(db_dir,outputdir):
	if not os.path.exists(db_dir):
		logger.error("Database does not exist at {}".format(db_dir))
		return
	orgs = os.walk(db_dir).next()[1] # get organism names
	xml_path = os.path.join(outputdir, "grsnp_gfs.xml")
	with open(xml_path,"wb") as writer:
		for o in orgs:
			blacklist = []
			# read in names of tracks to ignore
			blacklist_path = os.path.join(db_dir,o,"blacklist.txt")
			if os.path.exists(blacklist_path):
				with open(blacklist_path) as f:
					blacklist = [line.strip() for i,line in enumerate(f)]

			# generate the xml file for galaxy's checkbox tree
			writer.write("""<filter type="data_meta" data_ref="background" meta_key="dbkey" value="{}">\n\t<options>\n""".format(o))
			tmp = dir_as_xml(os.path.join(db_dir,o),blacklist).split("\n")
			tmp = "\n".join(tmp[1:-2]) # remove the first 'option' entry as this is the organism directory
			writer.write(tmp + "</options>\n</filter>") 
	logger.info("Created galaxy xml file {}".format(xml_path))		


def _exclude(x):
	ls = ["Raw","Align","Signal"]
	for l in ls:
		if l in x:
			return True
	if x.startswith("chr"):
		return True
	return False


def dir_as_xml(path, blacklist):
	''' Code adapted from
	 from: http://stackoverflow.com/questions/2104997/os-walk-python-xml-representation-of-a-directory-structure-recursion
	'''
	result = '<option name={} value={}>\n'.format(xml_quoteattr(os.path.basename(path))
													,xml_quoteattr(path))
	for item in os.listdir(path):
		itempath = os.path.join(path, item)
		if os.path.isdir(itempath):
			result += '\n'.join('  ' + line for line in 
			dir_as_xml(os.path.join(path, item),blacklist).split('\n'))
		elif os.path.isfile(itempath) and itempath.endswith(".bed.gz") and base_name not in blacklist:
			result += '  <option name={} value={}/>\n'.format(xml_quoteattr(base_name(item)), xml_quoteattr(os.path.join(path,item)))
	result += '</option>\n'
	return result

if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog="python -m grsnp.dbcreator", description='Creates the GenomeRunner SNP Database. Example: python -m grsnp.dbcreator -d /home/username/grsnp_db/ -g mm9', epilog='IMPORTANT: Execute DBCreator from the database folder, e.g., /home/username/grsnp_db/. Downloaded files from UCSC are placed in ./downloads database created in ./grsnp_db.')
	parser.add_argument("--data_dir" , "-d", nargs="?", help="Set the directory where the database to be created. Use absolute path. Example: /home/username/grsnp_db/. Required", required=True)
	parser.add_argument('--organism','-g', nargs="?", help="The UCSC code of the organism to use for the database creation. Default: hg19 (human). Required", default="hg19")
	parser.add_argument('--featurenames','-f', nargs="?", help='The names of the specific genomic feature tracks to create, comma separated (Example: knownGene, evoFold)', default="")
	parser.add_argument('--max','-m', nargs="?", help="Limit the number of features to be created within each group.",type=int)
	parser.add_argument('--galaxy', help="Create the xml files needed for Galaxy. Outputted to the current working directory.", action="store_true")
	parser.add_argument('--score', '-s', help="Commas separated list of score percentiles.", nargs='?',default="")
	parser.add_argument('--scoreonly','-o', help="Only filter by score. Skips downloading and installing new GFs.", action="store_true")



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
		
	global ftp, max_install_num
	ftp = ftplib.FTP(ftp_server, timeout=1800) # Connection timeout 0.5h
	ftp.login(username,password)
	outputdir=os.path.join(args["data_dir"],'grsnp_db')

	# if scoreonly is passed, then skip adding new GFs
	if not args['scoreonly']:
		if args['galaxy']:
			usrdir = raw_input("Enter directory of Galaxy containing the run.sh file. If left blank, grsnp_gfs.xml file will be outputted in the cwd: \n")
			if usrdir == '': 
				usrdir = os.getcwd()
			else:
				usrdir = os.path.join(usrdir,"tool-data")
				if not os.path.exists(usrdir):
					logger.error("Galaxy tool-data does not exist")
					sys.exist()
			create_galaxy_xml_files(outputdir,usrdir)		
			sys.exit()
		if args['organism'] is not None: # Only organism is specified. Download all organism-specific features
			gfs = args["featurenames"].split(",")
			trackdbpath = download_trackdb(args['organism'],outputdir)
			create_feature_set(trackdbpath,args['organism'],args["max"],gfs)			
		else:
			print "ERROR: Requires UCSC organism code.  Use --help for more information"
			sys.exit()

	# load score from minmax.txt file created earlier
	minmax = load_minmax(os.path.join(outputdir,args['organism'],"minmax.txt"))		

	### Second Step: Create subdirectories for score and filter data by score percentile
	# create sub directories for score percentiles and populate with score-filtered GF data
	# gather all directories (groups) in the database
	print "Filtering GFs by Score"
	orgdir = os.path.join(outputdir,args['organism'])
	dirs = [name for name in os.listdir(orgdir)
		if os.path.isdir(os.path.join(orgdir, name))]
	for d in dirs:
		# gather all paths
		gfs = []
		for base, tmp, files in os.walk(os.path.join(orgdir,d)):
				gfs += [os.path.join(base,f) for f 	 
					in files if f.endswith(('.gz', '.bb'))]
		for gf_path in gfs: 
			for pct_score in args['score']:		
				[score_min,score_max] = minmax[base_name(gf_path)].split(",")
				# calculate threshold score
				if score_min == 'NA':
					continue
				score_min,score_max = float(score_min),float(score_max) 
				thresh_score = score_min + (score_max-score_min)*float(pct_score)/100
				logger.info("MinMax stats for {}: Min={}, Max={}, {} pct_thresh={}".format(base_name(gf_path), score_min,score_max,pct_score,thresh_score))
				# is this safe? It searches /dirpath/grsnp_db/subdirs/gf.txt and replaces /grsnp_db/ with /grsnp_db_[score]/
				gf_path_out =gf_path.replace('/grsnp_db/','/grsnp_db_{}/'.format(pct_score))
				if not os.path.exists(os.path.split(gf_path_out)[0]):
					os.makedirs(os.path.split(gf_path_out)[0])
				gf_path_out_woext = os.path.join(os.path.split(gf_path_out)[0],base_name(gf_path_out))
				filter_by_score(gf_path, gf_path_out_woext,thresh_score)


	root_dir = os.path.dirname(os.path.realpath(__file__))
	readme = open(os.path.join(root_dir,"grsnp_db_readme.txt")).read()
	with open("grsnp_db_readme.txt","wb") as writer:
		writer.write(readme)
	print "FINISHED: Downloaded files from UCSC are placed in {}.  Database created in {}".format(os.path.join(args["data_dir"],"downloads"),os.path.join(args["data_dir"],"grsnp_db`"))


