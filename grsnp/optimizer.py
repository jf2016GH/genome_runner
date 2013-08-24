import sys
import logging
from logging import FileHandler,StreamHandler
import os
import gzip
import subprocess
import argparse
import hypergeom4 as hpgm
import pdb

# connection information for the ucsc ftp server
logger = logging.getLogger()




def create_bkg_gf_overlap_db(gf_dir,background_dir):
	""" Used to precalculate the overlapStatistics for the GFs against each of the
	default backgrounds.
	"""

	gf_bg_stats,list_completed_gfs = {},[]
	all_gfs = []
	backgrounds,gfs=[],[]
	db_path = os.path.join(gf_dir,"bkg_overlaps.gr")
	db_path_tmp = db_path + ".tmp"


	# Read in all completed GF in the partially completed bkg_overlaps file
	if os.path.exists(db_path_tmp):
		list_completed_gfs = [x.split("\t")[0] for x in open(db_path_tmp).read().split("\n")] 

	# gather all directories (groups) in the database
	dirs = [name for name in os.listdir(gf_dir)
		if os.path.isdir(os.path.join(gf_dir, name))]

	# Gather backgrounds paths
	backgrounds = [os.path.join(background_dir, f) for f in os.listdir(background_dir) if f.endswith(('.gz', '.bb',".txt",".bed"))]

	cur_prog,prog_max = 1,_count_gfs(gf_dir)
	# Process each category of GFs
	for d in dirs:
		logger.info("Running overlapStatistics for all GFs in {}".format(d))
		# Gather gfs paths
		gfs = []
		for base, d, files in os.walk(os.path.join(gf_dir,d)):
				gfs += [os.path.join(base, f) for f 
					in files if f.endswith(('.gz', '.bb'))]
		# create empty list entries for each gf
		for g in gfs:
			gf_bg_stats[g] = []
		all_gfs += gfs
		prog_gf  = 1
		# Run overlap analysis for each GF (g) and save to results in dictionary with key equal to 'g'
		for g in gfs:
			if g not in list_completed_gfs:	
				_write(g+"\t",db_path_tmp)
				for bg in backgrounds:
					res = hpgm.get_overlap_statistics(g,[bg])
					_write(os.path.join(background_dir,res[0]["queryfile"])+":"+str(res[0]["intersectregions"])+":"+str(res[0]["queryregions"])+",",db_path_tmp)  # write out bg_obs
				logger.info("Processed {}".format(g))
				print "Processed {} of {}.".format(cur_prog,prog_max)
				_write("\n",db_path_tmp)
			else:
				logger.info("Stats already exist for {}, skipping.".format(g))
			cur_prog += 1

	if os.path.exists(db_path): os.remove(db_path) # remove old bgs_overlap
	os.rename(db_path_tmp,db_path) # release new bgs_overlap to GRSNP
	logger.info("Completed.")


def _count_gfs(grsnp_db):
	x = 0
	for root, dirs, files in os.walk(grsnp_db):
		for f in files:
			if f.endswith(".bed.gz"):
				x = x+1
	return x

def _write(line,path):
	with open(path,"a") as wr:
		wr.write(line)	


if __name__ == "__main__":
	global logger
	parser = argparse.ArgumentParser(description="""Pre calculates the overlapStatistics for each of the backgrounds in /custom_data/backgrounds/[organism] 
												 and genomic features in /grsnp_db/[organism]. Creates a file in the database called called bkg_overlap.gr which is used
												 by the server in hypergeom4.""")
	parser.add_argument('--organism','-g',required=True, help='The UCSC code of the organism to be downloaded (example: hg19 (human))')
	parser.add_argument('--data_dir','-d', required=True, help="Directory containing of the GRSNP database. Use ABSOLUTE path.")
	args = vars(parser.parse_args())


	hdlr = logging.FileHandler(os.path.join(args["data_dir"],'genomerunner_dbcreator.log'))
	hdlr_std = StreamHandler()
	formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
	hdlr.setFormatter(formatter)
	hdlr_std.setFormatter(formatter)
	logger.addHandler(hdlr)
	logger.addHandler(hdlr_std)
	logger.setLevel(logging.INFO)

	# Ask if use wants to continue partially run optimization
	path_tmp = os.path.join(args["data_dir"],"grsnp_db",args["organism"],"bkg_overlaps.gr") + ".tmp"
	if os.path.exists(path_tmp):
		in_var = raw_input("Temporary file exists at {}. Do you want to continue (yes) or start from scratch (no)?".format(path_tmp))
		if in_var.lower() == "no": 
			os.remove(path_tmp)

	background_dir = os.path.join(args["data_dir"],"custom_data","backgrounds",args["organism"])
	gfs_dir = os.path.join(args["data_dir"],"grsnp_db",args["organism"])
	if not os.path.exists(gfs_dir):
		print "ERROR: grsnp_db does not exist.  Use grsnp.dbcreator to create a database."
		sys.exit()
	if not os.path.exists(background_dir):
		print "ERROR: No backgrounds found in default background directory {}.  Please add backgrounds.".format(background_dir)
		sys.exit()
	logger.info("Pre calculating statistics for GR database in")
	create_bkg_gf_overlap_db(gf_dir=gfs_dir,background_dir=background_dir)

