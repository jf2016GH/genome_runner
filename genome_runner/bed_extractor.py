import sys
import os
import ftplib
import sqlite3
import string
from contextlib import closing
import subprocess
import argparse

# connection information for the uscs ftp server
server = 'hgdownload.cse.ucsc.edu'
directory = '/goldenPath/{}/database'
username = 'anonymous'
password = ''
downloadtodir = ''

def _get_columns()
''' Extracts the column names, in order, from the provided sql file
	Variables must be line seperated for the function to work
'''
def create_database(organism):
	if os.path.exists("data") == False:
		try:
			os.makedirs("data")
		except:
			pass
	dbpath = os.path.join("data",organism + ".db3")

	# download the trackdb.sql file
	track_sql_path = download_trackdb(organism)	
	print "creating database and adding trackdb"
	args = ["sqlite3", organism + ".db3",_convert_mysqlquery_sqlite(track_sql_path)]
	subprocess.call(args)
	print "created database at " 


		
def execute(dbpath, qry):
	with closing(sqlite3.connect(dbpath)) as conn:
		with closing(conn.cursor()) as c:
			c.execute(qry)
			return [r for r in c]

