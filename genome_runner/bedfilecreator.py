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
		sw.write(text.replace('\\\n','<br />').replace('\\\t','     '))
	return sqloutputpath
	



def extract_bed6(outputpath,datapath,colnames):
	colstoextract = ['chrom','chromStart','chromEnd','name','score','strand']
	print "Outpath is: {}".format(outputpath)
	with gzip.open(datapath) as dr:
		with gzip.open(outputpath,"wb") as bed:
			while True:
				line = dr.readline()
				if line == "":
					break
				r  = dict(zip(colnames,line.split('\t')))
				row = []
				for col in colnames:
					row.append(r[col]) 
				bed.write("\t".join(map(str,row))+"\n")


def extract_bed3(outputpath,datapath,colnames):
	colstoextract = ['chrom','chromStart','chromEnd']
	with gzip.open(outputpath,"wb") as bed:
		with gzip.open(datapath) as dr:
			while True:
				line = dr.readline()
				if line == "":
					break
				r  = dict(zip(colnames,line.split('\t')))
				row = [r["chrom"],r["chromStart"],r["chromEnd"],".",".","."]
				bed.write("\t".join(map(str,row))+"\n")

def extract_genepred(outputpath,datapath,colnames):
	colstoextract = ['chrom','txStart','txEnd','name','strand']
	exonpath = outputpath.split(".")[0]+"_exon.gz"
	with gzip.open(datapath) as dr:
		with gzip.open(outputpath,"wb") as bed:
			with gzip.open(exonpath+".temp","wb") as exonbed:
				while True:
					line = dr.readline()
					if line == "":
						break
					r = dict(zip(colnames,line.split('\t')))
					# extract the gene data inserts a blank for score
					row = [r['chrom'],r['txStart'],r['txEnd'],r['name'],'.',r['strand']]
					bed.write("\t".join(map(str,row))+"\n")
					# extract the exon data
					for (s,e) in zip(r["exonStarts"].split(","),r["exonStarts"].split(",")):
						if s != '':
							rowexon = [r['chrom'],s,e,r['name'],'.',r['strand']]
							exonbed.write("\t".join(map(str,rowexon))+"\n")
	# remove the .temp extension from the exon file 
	os.rename(exonpath+".temp",exonpath)

	

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
			if numdownloaded[row["type"]] <=5000000:
				sqlpath = download_uscs_file(organism,row["tableName"] + ".sql","downloads")
				download_uscs_file(organism,row["tableName"] + ".txt.gz","downloads")
				if sqlpath != '':
					print "converting",row['tableName'], " into proper bed format"
					try:
						outpath = os.path.join(outputdir,row["grp"],'Tier' + row["visibility"],row["tableName"]+".gz")
						if not os.path.exists(os.path.dirname(outpath)):
							os.makedirs(os.path.dirname(outpath))
						if os.path.exists(outpath) == False:
							# removes the .temp file, to prevent duplicate data from being written
							if os.path.exists(outpath+".temp"):
								os.remove(outpath+".temp")
							# converts the uscs data into propery bed format
							print "Converting into proper bed format. {}".format(os.path.splitext(sqlpath)[0] + ".txt.gz")
							preparebed[row["type"]](outpath+".temp",os.path.splitext(sqlpath)[0]+".txt.gz",get_column_names(os.path.splitext(sqlpath)[0]+".sql"))
							# remove the .temp file extension to activate the GF
							os.rename(outpath+".temp",outpath)
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

