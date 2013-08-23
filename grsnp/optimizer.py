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

	gf_bg_stats = {}
	all_gfs = []
	backgrounds,gfs=[],[]
	dirs = [name for name in os.listdir(gf_dir)
		if os.path.isdir(os.path.join(gf_dir, name))]

	print "DRIS: ",dirs

	# Gather backgrounds paths
	backgrounds = [os.path.join(background_dir, f) for f in os.listdir(background_dir) if f.endswith(('.gz', '.bb',".txt",".bed"))]

	# Process each category of GFs
	for d in dirs:
		logger.info("Running overlapStatistics for all GFs in {}".format(d))
		# Gather gfs paths
		gfs = []
		for base, dirs, files in os.walk(os.path.join(gf_dir,d)):
				gfs += [os.path.join(base, f) for f 
					in files if f.endswith(('.gz', '.bb'))]
		# create empty list entries for each gf
		for g in gfs:
			gf_bg_stats[g] = []
		all_gfs += gfs
		# Run overlap analysis for each GF (g) and save to results in dictionary with key equal to 'g'
		for g in gfs:	
			for bg in backgrounds:
				print g	
				gf_bg_stats[g] += hpgm.get_overlap_statistics(g,[bg])
	# output overlap statistics
	dp_path = os.path.join(gf_dir,"bkg_overlaps.gr")
	logger.info("Outputing results to {}".format(dp_path))
	with open(dp_path,"wb") as out:
		for g in all_gfs: # cycle through all of the GFs in the database
			out.write(g+"\t")
			if g in gf_bg_stats:	# check if the current GF was run successfully	
				for res in gf_bg_stats[g]: # lookup the results for the current GF
					out.write(os.path.join(background_dir,res["queryfile"])+":"+str(res["intersectregions"])+":"+str(res["queryregions"])+",")
				out.write("\n")


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
	logger.addHandler(hdlr)
	logger.addHandler(hdlr_std)
	logger.setLevel(logging.INFO)

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

