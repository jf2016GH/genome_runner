import sys
import logging
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
# connection information for the uscs ftp server
server = 'hgdownload.cse.ucsc.edu'
directory = '/goldenPath/{}/database'
username = 'anonymous'
password = ''
eventlog = logging.getLogger("grdbcreator")
logfilename = "gr_log.txt"
			
# downloads the specified file from uscs.  Saves it with a .temp extension untill the download is complete.
def download_uscs_file(organism,filename,downloaddir):
	''' Downloads the filename from the USCS ftp server and saves it
	in a folder with the same name as the organism.
	'''
	outputpath = ''
	if downloadtodir != None and downloadtodir != '':
		outputdir = os.path.join(downloaddir,organism)
	else:
		outputdir = organism
	try:
		if os.path.exists(outputdir) == False and outputdir != '':
			print "creating directory {}".format(outputdir)
			os.makedirs(outputdir)
	except Exception, e:
		print e
		eventlog.warning("ERROR: could not create folder at {} for {}".format(outputdir,filename),logfilename)
		return '' 
	
	try:
		outputpath = os.path.join(outputdir,filename)
		if not os.path.exists(outputpath):
			with open(outputpath + ".temp",'wb') as fhandle:  
				print 'Downloading {} from USCS'.format(filename)
				ftp = ftplib.FTP(server)
				ftp.login(username,password)
				ftp.cwd(directory.format(organism))
				ftp.retrbinary('RETR ' + "{}".format(filename),fhandle.write)
				os.rename(outputpath+".temp",outputpath)
				print 'Finished downloading {} from USCS'.format(filename)
		else:
			print '{} already exists, skipping download'.format(outputpath)
	except Exception, e:
		print e
		eventlog.warning("ERROR: could not download the {} sql file. Names ARE case sensitive.".format(filename),logfilename)
		return '' 

	return outputpath 


def download_trackdb(organism):
	''' Downloads the trackdb.sql and trackDb.txt.gz from the USCS ftp server and saves it in a folder with the same name as the organism.
		Returns the path of the downloaded .sql file
	'''
	sqloutputpath = download_uscs_file(organism,"trackDb.sql",downloadtodir)
	dataoutpath = download_uscs_file(organism,"trackDb.txt.gz",downloadtodir)
	' replace all of the \\\n characters in the html column with <br />'
	text = gzip.open(dataoutpath).read()
	with gzip.open(dataoutpath,'wb') as sw:
		sw.write(text.replace('\\\n','<br />'))

	return sqloutputpath
	



def extract_bed6(outputpath,tabledata):
	colstoextract = ['chrom','chromStart','chromEnd','name','score','strand']
	data = []
	for r in tabledata:
		row = []
		for col in colstoextract:
			row.append(r[col])
		data.append(row)

	write_bed(outputpath,data)

def extract_bed3(outputpath,tabledata):
	colstoextract = ['chrom','chromStart','chromEnd']
	data = []
	for r in tabledata:
		row = []
		for col in colstoextract:
			row.append(r[col])
		data.append(row)

	write_bed(outputpath,data)

def extract_genepred(outputpath,tabledata):
	colstoextract = ['chrom','txStart','txEnd','name','strand']
	genedata = []
	exondata = []
	for r in tabledata:
		# extract the gene data inserts a blank for score
		genedata.append((r['chrom'],r['txStart'],r['txEnd'],r['name'],'.',r['strand']))
		# extract the exon data
		for (s,e) in zip(r["exonStarts"].split(","),r["exonStarts"].split(",")):
			if s != '':
				exondata.append((r['chrom'],s,e,r['name'],'.',r['strand']))
	
	write_bed(outputpath,genedata)
	# write the exon data
	write_bed(outputpath.split(".")[0] + "_exon.gz",exondata)

# writes the formated bed file to a .gz.temp.  After the write is complete, it is renamed to .gz
def write_bed(outputpathgz,data):
	''' writes the data to a .gz bed data file
	'''
	print 'writting the data to {}'.format(outputpathgz)
	if os.path.exists(os.path.dirname(outputpathgz)) != True:
		os.makedirs(os.path.dirname(outputpathgz))	
	
	with gzip.open(outputpathgz + ".temp",'wb') as bed:
		for row in data:
			bed.write("\t".join(map(str,row))+"\n")
	os.rename(outputpathgz+".temp",outputpathgz)

def get_column_names(sqlfilepath):
	''' extracts the column names from the .sql file and returns them
	'''
	tdbsql = open(sqlfilepath).read()
	# creates a tuple containing the column names
	tbdcolumns = re.findall("\n\s+`(.+?)`", tdbsql, re.DOTALL)
	return tbdcolumns	
preparebed = {"bed 6" : extract_bed6, "bed 6 +" : extract_bed6, "bed 3" : extract_bed3,"genePred" : extract_genepred}
numdownloaded = {"bed 6" : 0, "bed 6 +" : 0, "bed 3" : 0, "genePred" : 0}


def download_bedfiles(trackdbpath,organism):
	outputdir = os.path.dirname(trackdbpath)
	trackdb = load_tabledata_dumpfiles(os.path.splitext(trackdbpath)[0])
	prog, num = 0,len(trackdb) 
	for row in trackdb:	
		print 'Processing files {} of {}'.format(prog,num)
		if row['type'] in preparebed:
			# DEBUG this line limits the number of GRF to download
			for d in numdownloaded.values():
				print d
			if numdownloaded[row["type"]] <= 10:
				sqlpath = download_uscs_file(organism,row["tableName"] + ".sql","downloads")
				download_uscs_file(organism,row["tableName"] + ".txt.gz","downloads")
				if sqlpath != '':
					data = load_tabledata_dumpfiles(os.path.splitext(sqlpath)[0])
					print "converting",row['tableName'], " into proper bed format"
					try:
						outpath = os.path.join(outputdir,row["grp"],'Tier' + row["visibility"],row["tableName"]+".gz")
						if os.path.exists(outpath) == False:
							preparebed[row["type"]](outpath,data)
						else:
							print "{} already exists, skipping extraction".format(outpath)
						numdownloaded[row["type"]] += 1
					except Exception, e:
						print "ERROR unable to convert {} into bed".format(row["tableName"])
						print e
						continue
		prog += 1
	return "created database"
	



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
	parser = argparse.ArgumentParser(description='Create the GenomeRunner Database')
	parser.add_argument('--output','-o', help='The folder to download the data files to')
	parser.add_argument('--organism','-g', help='The USCS code of the organism to be downloaded (example: hg19 (human))')
	parser.add_argument('--filename','-f', help='The name of the specific file to download (example: trackDb.sql). If provided, only this file will be downloaded')
	args = vars(parser.parse_args())
	downloadtodir='download'
	if args['output'] is not None:
		downloadtodir = args['output']
	if args['organism'] is not None and args['filename'] is not None:
		download_uscs_file(args['organism'],args['filename'],downloadtodir)
	elif args['organism'] is not None and args['filename'] is None:
		trackdbpath = download_trackdb(args['organism'])
		download_bedfiles(trackdbpath,args['organism'])
	else:
		trackdbpath = download_trackdb('hg19')
		print trackdbpath
		download_bedfiles(trackdbpath,'hg19')

