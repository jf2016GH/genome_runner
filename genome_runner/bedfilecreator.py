import sys
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

# connection information for the ucsc ftp server
server = 'hgdownload.cse.ucsc.edu'
directory = '/goldenPath/{}/database'
username = 'anonymous'
password = ''

logger = logging.getLogger('genomerunner.dbcreator')
hdlr = logging.FileHandler('genomerunner_dbcreator.log')
hdlr_std = StreamHandler()
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.addHandler(hdlr_std)
logger.setLevel(logging.INFO)

ftp = ftplib.FTP(server)
ftp.login(username,password)
			
# downloads the specified file from ucsc.  Saves it with a .temp extension untill the download is complete.
def download_ucsc_file(organism,filename,downloaddir):
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
				FTP.quit()
				os.rename(outputpath+".temp",outputpath)
				logger.info( 'Finished downloading {} from UCSC'.format(filename))
		else:
			logger.info( '{} already exists, skipping download'.format(outputpath))
	except Exception, e:
		logger.warning( e)
		logger.warning("Could not download the {} sql file. Names ARE case sensitive.".format(filename))
		return '' 

	return outputpath 


def download_trackdb(organism,outputdir):
	''' Downloads the trackdb.sql and trackDb.txt.gz from the UCSC ftp server and saves it in a folder with the same name as the organism.
		Returns the path of the downloaded .sql file
	'''
	sqloutputpath = download_ucsc_file(organism,"trackDb.sql",outputdir)
	dataoutpath = download_ucsc_file(organism,"trackDb.txt.gz",outputdir)
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
	colstoextract = ['chrom','chromStart','chromEnd','name','score','strand']
	# Checks if all of the columns exist in the table.  If not extract_bed5 is tried instead
	if _check_cols(colnames,colstoextract):

		logger.info( "Outpath is: {}".format(outputpath))
		with gzip.open(datapath) as dr:
			with gzip.open(outputpath,"wb") as bed:
				while True:
						line = dr.readline().strip('\r').rstrip('\n')
						if line == "":
							break
						r  = dict(zip(colnames,line.split('\t')))
						row = []
						row = [r["chrom"],r["chromStart"],r["chromEnd"],r["name"],r["score"] if r["score"] != "." else "0",r["strand"]]
						bed.write("\t".join(map(str,row))+"\n")
	else:
		logger.warning("Nonstandard bed6, attempting extraction as bed5")
		extract_bed5(outputpath,datapath,colnames)

def extract_bed5(outputpath,datapath,colnames):
	colstoextract = ['chrom','chromStart','chromEnd','name','score']
	if _check_cols(colnames,colstoextract):
		with gzip.open(outputpath,"wb") as bed:
			with gzip.open(datapath) as dr:
				while True:
					line = dr.readline().rstrip('\r').rstrip('\n')
					if line == "":
						break
					r  = dict(zip(colnames,line.split('\t')))
					row = [r["chrom"],r["chromStart"],r["chromEnd"],r["name"],r["score"] if r["score"] != "." else "0","."]
					bed.write("\t".join(map(str,row))+"\n")
	else:
		logger.warning("Nonstandard bed5, attempting extraction as bed4")
		extract_bed4(outputpath,datapath,colnames)
	
def extract_bed4(outputpath,datapath,colnames):
	colstoextract = ['chrom','chromStart','chromEnd','name']
	if _check_cols(colnames,colstoextract):
		with gzip.open(outputpath,"wb") as bed:
			with gzip.open(datapath) as dr:
				while True:
					line = dr.readline().rstrip('\r').rstrip('\n')
					if line == "":
						break
					r  = dict(zip(colnames,line.split('\t')))
					row = [r["chrom"],r["chromStart"],r["chromEnd"],r["name"],"0","."]
					bed.write("\t".join(map(str,row))+"\n")
	else:
		logger.warning("Nonstandard bed4, attempting extraction as bed3")
		extract_bed3(outputpath,datapath,colnames)

def extract_bed3(outputpath,datapath,colnames):
	colstoextract = ['chrom','chromStart','chromEnd']
	with gzip.open(outputpath,"wb") as bed:
		with gzip.open(datapath) as dr:
			while True:
				line = dr.readline().strip('\r').rstrip('\n')
				if line == "":
					break
				r  = dict(zip(colnames,line.split('\t')))
				row = [r["chrom"],r["chromStart"],r["chromEnd"],".","0","."]
				bed.write("\t".join(map(str,row))+"\n")

def extract_genepred(outputpath,datapath,colnames):
	colstoextract = ['chrom','txStart','txEnd','name','strand']
	exonpath = outputpath.split(".")[0]+"_exon.gz"
	with gzip.open(datapath) as dr:
		with gzip.open(outputpath,"wb") as bed:
			from os.path import basename
			with gzip.open(exonpath+".temp","wb") as exonbed:
				while True:
					line = dr.readline().rstrip('\r').rstrip('\n')
					if line == "":
						break
					r = dict(zip(colnames,line.split('\t')))
					# extract the gene data inserts a blank for score
					row = [r['chrom'],r['txStart'],r['txEnd'],r['name'],'0',r['strand']]
					bed.write("\t".join(map(str,row))+"\n")
					# extract the exon data
					for (s,e) in zip(r["exonStarts"].split(","),r["exonEnds"].split(",")):
						if s != '':
							rowexon = [r['chrom'],s,e,r['name'],'0',r['strand']]
							exonbed.write("\t".join(map(str,rowexon))+"\n")
	# remove the .temp extension from the exon file 
	os.rename(exonpath+".temp",exonpath)

def extract_rmsk(outputpath,datapath,colnames):
	colstoextract = ['genoName','genoStart','genoEnd','repClass', 'strand','swScore']
	with gzip.open(outputpath,"wb") as bed:
		with gzip.open(datapath) as dr:
			while True:
				line = dr.readline().strip('\r').rstrip('\n')
				if line == "":
					break
				r = dict(zip(colnames,line.split('\t')))
				row = [r["genoName"],r["genoStart"],r["genoEnd"],r["repClass"],r["strand"],r["swScore"]]
				bed.write("\t".join(map(str,row))+"\n")

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
				"broadPeak": extract_bed6,
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
				"rmsk" : extract_rmsk,
				"factorSource" : extract_bed6}
				
numdownloaded = collections.defaultdict(int)

def create_feature_set(trackdbpath,organism):
	outputdir = os.path.dirname(trackdbpath)
	trackdb = load_tabledata_dumpfiles(os.path.splitext(trackdbpath)[0])
	prog, num = 0,len(trackdb)
	added_features = [] 
	notsuptypes, outpath = set([]),""
	for row in trackdb:
		logger.info( 'Processing files {} of {}'.format(prog,num))
		if row['type'] in preparebed:
			# DEBUG this line limits the number of GRF to download
			if numdownloaded[row["type"]] <= 5:
				sqlpath = download_ucsc_file(organism,row["tableName"] + ".sql","downloads")
				download_ucsc_file(organism,row["tableName"] + ".txt.gz","downloads")
				if sqlpath != '':
					logger.info( "converting"+row['tableName']+ " into proper bed format")
					try:
						outpath = os.path.join(outputdir,row["grp"],'Tier' + row["visibility"],row["tableName"]+".gz")
						if not os.path.exists(os.path.dirname(outpath)):
							os.makedirs(os.path.dirname(outpath))
						if os.path.exists(outpath) == False:
							# removes the .temp file, to prevent duplicate data from being written
							if os.path.exists(outpath+".temp"):
								os.remove(outpath+".temp")
							# converts the ucsc data into propery bed format
							logger.info( "Converting into proper bed format. {}".format(os.path.splitext(sqlpath)[0] + ".txt.gz"))
							preparebed[row["type"]](outpath+".temp",os.path.splitext(sqlpath)[0]+".txt.gz",get_column_names(os.path.splitext(sqlpath)[0]+".sql"))
							# remove the .temp file extension to activate the GF
							os.rename(outpath+".temp",outpath)
							added_features.append(outpath)
						else:
							logger.info( "{} already exists, skipping extraction".format(outpath))							
						numdownloaded[row["type"]] += 1
					except Exception, e:
						exc = trace.format_exc()
						logger.warning( "Unable to convert {} into bed".format(row["tableName"]))
						logger.warning(exc)
						continue
		else:
			if 'big' not in row['type']:
				notsuptypes.add(row['type'])

		prog += 1
		# cleanup the temporary files
		if os.path.exists(outpath + ".temp"): os.remove(outpath+".temp")

	logger.info( "The following types are not supported (includes all 'big' file types):\n " + str(notsuptypes))
	logger.info("The following features were added to the database: \n{}".format(added_features))
	logger.info("A count of features added by type: ")
	for k,d in numdownloaded.iteritems():
		logger.info( k + ":" + str(d))
	return "created database"
	
def create_single_feature(trackdbpath,organism,feature):
	''' Downloads a single feature and adds it to the genomerunner flat file database'''

	outputdir = os.path.dirname(trackdbpath)
	trackdb = load_tabledata_dumpfiles(os.path.splitext(trackdbpath)[0])
	f_type = _gettype(feature,trackdb)
	f_info = _get_info(feature,trackdb)
	# is the feature in trackDb
	if f_info != False:
		# is the feature type supported by the dbcreator
		if  f_type in preparebed:
			sqlpath = download_ucsc_file(organism,f_info["tableName"] + ".sql","downloads")
			download_ucsc_file(organism,f_info["tableName"] + ".txt.gz","downloads")
			if sqlpath != '':
				logger.info( "converting"+f_info['tableName']+ " into proper bed format")
				try:
					outpath = os.path.join(outputdir,f_info["grp"],'Tier' + f_info["visibility"],f_info["tableName"]+".gz")
					if not os.path.exists(os.path.dirname(outpath)):
						os.makedirs(os.path.dirname(outpath))
					# if the feature is not in the database, add it
					if os.path.exists(outpath) == False:
						# removes the .temp file, to prevent duplicate data from being written
						if os.path.exists(outpath+".temp"):
							os.remove(outpath+".temp")
						# converts the ucsc data into propery bed format
						logger.info( "Converting into proper bed format. {}".format(os.path.splitext(sqlpath)[0] + ".txt.gz"))
						preparebed[f_info["type"]](outpath+".temp",os.path.splitext(sqlpath)[0]+".txt.gz",get_column_names(os.path.splitext(sqlpath)[0]+".sql"))
						# remove the .temp file extension to activate the GF
						os.rename(outpath+".temp",outpath)
					else:
						logger.info( "{} already exists, skipping extraction".format(outpath))
					numdownloaded[f_info["type"]] += 1
				except Exception, e:
					exc = trace.format_exc()
					logger.warning( "Unable to convert {} into bed".format(f_info["tableName"]))
					logger.warning(exc)
		else:
			logger.warning("{} is a type {}, which is not supported".format(feature,f_type))
	else:
		logger.warning( "Could not find {} in trackDb".format(feature))

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


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Creates the GenomeRunner Database.  Downloaded files from UCSC are placed in /Downloads.  Converted files are placed in /Released')
	parser.add_argument('--organism','-g', help='The UCSC code of the organism to be downloaded (example: hg19 (human))')
	parser.add_argument('--featurename','-f', help='The name of the specific genomic feature track to create (example: knownGene)')
	args = vars(parser.parse_args())
	outputdir='released'
	if args['organism'] is not None and args['featurename'] is None: # Only organism is specified. Download all organism-specific features
		trackdbpath = download_trackdb(args['organism'],outputdir)
		#pdb.set_trace()
		create_feature_set(trackdbpath,args['organism'])
	elif args['organism'] is not None and args['featurename'] is not None: # Both organism and feature name are specified. Download this feature for a given organism
		trackdbpath = download_trackdb(args['organism'],outputdir)
		create_single_feature(trackdbpath,args['organism'],args['featurename'])
	elif args['organism'] is None and args['featurename'] is not None: # Warning in case of only feature name is supplied
		print "To add a specific feature to the local database, please supply an organism assembly name"
	else:
		print "Requires UCSC organism code.  Use --help for more information"
