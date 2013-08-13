import sys
import logging
from logging import FileHandler,StreamHandler
import os
import gzip
import subprocess
import argparse

# connection information for the ucsc ftp server
logger = logging.getLogger('genomerunner.dbcreator')
hdlr = logging.FileHandler('genomerunner_dbcreator.log')
hdlr_std = StreamHandler()
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.addHandler(hdlr_std)
logger.setLevel(logging.INFO)





def create_bkg_gf_overlap_db(organism, gf_dir,background_dir):
	""" Used to precalculate the overlapStatistics for the GFs against each of the
	default backgrounds.
	"""

	background_stats = {}
	all_gfs = []
	backgrounds,gfs=[],[]
	gf_dir = os.path.join(gf_dir,organism)
	background_dir = os.path.join(background_dir,organism)
	dirs = [name for name in os.listdir(gf_dir)
		if os.path.isdir(os.path.join(gf_dir, name))]

	print "DRIS: ",dirs

	# Gather backgrounds paths
	backgrounds = [os.path.join(background_dir, f) for f in os.listdir(background_dir) if f.endswith(('.gz', '.bb',".txt",".bed"))]
	for b in backgrounds:
		background_stats[b] = []

	# Process each category of GFs
	for d in dirs:
		logger.info("Running overlapStatistics for all GFs in {}".format(d))
		# Gather gfs paths
		gfs = []
		for base, dirs, files in os.walk(os.path.join(gf_dir,d)):
				gfs += [os.path.join(base, f) for f 
					in files if f.endswith(('.gz', '.bb'))]
		print
		all_gfs += gfs
		# Run overlap analysis for each GF and save to results
		for b in backgrounds:		
			background_stats[b] += get_overlap_statistics(b,gfs)
		print background_stats
	# output overlap statistics
	dp_path = os.path.join(gf_dir,"bkg_overlaps.gr")
	logger.info("Outputing results to {}".format(dp_path))
	with open(dp_path,"wb") as out:
		for g in all_gfs: # output results for each GF in order
			out.write(g+"\t")		
			for k,v in background_stats.items():
				print v
				print "GF",g
				bk_obs = [x["intersectregions"] for x in v if x["queryfile"] == g]
				if len(bk_obs) != 0: out.write(k+":"+str(bk_obs[0])+",")
			out.write("\n")

def get_overlap_statistics(gf,fois):
    """Returns a dictionary with indicating how many hits exist for each foi against the gf
    gf: filepath for GF
    fois: list of FOI filepaths
    """
    results = []
    out = ""
    print "Background: ", [gf]
    print "GFs: ", fois
    out = subprocess.Popen(["overlapStatistics"] + [gf] + fois,stdout=subprocess.PIPE)
    out.wait()
    tmp = out.stdout.read()
    for x in tmp.split("\n")[1:]:
        if x != "":
            tmp = x.split("\t")
            foi_name,n,hit_count = tmp[0],tmp[2],tmp[3]
            results.append({"queryfile": foi_name,"queryregions": int(n),"intersectregions": int(hit_count)})
    return results


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='GenomeRunner GF database utilities')
	parser.add_argument('--organism','-g',required=True, help='The UCSC code of the organism to be downloaded (example: hg19 (human))')
	#parser.add_argument('--gf_dir','-d', required=True, help="Directory containing the files for all organisms. (i.e ./data for ./data/hg19)", action='store_true')
	#parser.add_argument('--background_dir','-b', required=True, help="Directory containing the backgrounds. (i.e ./backgrounds for ./backgrounds/hg19", action='store_true')

	args = vars(parser.parse_args())
	create_bkg_gf_overlap_db(args["organism"], gf_dir="released",background_dir="../backgrounds")

