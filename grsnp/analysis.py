#!/usr/bin/env python2
from __future__ import division

import argparse
import os



def main():
	global detailed_outpath, progress_outpath, run_files_dir, console_output, print_progress
	print_progress = True
	parser = argparse.ArgumentParser(
		description="Enrichment analysis of several sets of SNPs (FOIs) files against several genomic features (GFs). Example: python analysis.py foi_full_names.txt gf_full_names.txt /path_to_background/snp137.bed.gz")
	parser.add_argument("fois", nargs=1, help="Text file with paths to FOI files (unless -p used). Required")
	parser.add_argument("gfs", nargs=1,
						help="Text file with pathrs to GF files (unless -p used). GF files may be gzipped. Required")
	parser.add_argument("bg_path", nargs=1, help="Path to background, or population of all SNPs. Required")
	parser.add_argument("--run_annotation", "-a", help="Run annotation analysis", action="store_true")
	parser.add_argument("--run_files_dir", "-r", nargs="?",
						help="Set the directory where the results should be saved. Use absolute path. Example: /home/username/run_files/.",
						default=os.getcwd())
	parser.add_argument("--pass_paths", "-p",
						help="Pass fois and gfs as comma separated paths. Paths are saved in .fois and .gfs file.",
						action="store_true")
	parser.add_argument("--data_dir", "-d", nargs="?", type=str,
						help="Set the directory containing the database. Required for rsID conversion. Use absolute path. Example: /home/username/db_#.##_#.##.####/.",
						default="")
	parser.add_argument('--organism', '-g', nargs="?",
						help="The UCSC code of the organism to use. Required for rsID conversion. Default: hg19 (human).",
						default="hg19")
	default_test = "chisquare"
	parser.add_argument('--stat_test', '-s', nargs="?",
						help="Select the statistical test to use for calculating P-values. Default: {}. Available: chisquare, binomial, montecarlo_[# of simulations]".format(
							default_test), default=default_test)
	args = vars(parser.parse_args())
	if args['organism'] is None:
		print "--organism cannot be blank"
		return None
	if args['run_files_dir'] is None:
		print "--run_files_dir cannot be blank"
		return None
	if args["pass_paths"]:
		gf = args["gfs"][0].split(",")
		foi = args["fois"][0].split(",")
		if not os.path.exists(args["run_files_dir"]) and args['run_files_dir'] != "":
			os.mkdir(args['run_files_dir'])
		# write out the passed gf and foi paths into .gfs and .fois files.
		args["gfs"][0], args["fois"][0] = os.path.join(args["run_files_dir"], ".gfs"), os.path.join(
			args["run_files_dir"], ".fois")
		with open(args["gfs"][0], 'wb') as writer:
			writer.write("\n".join(gf))
		with open(args["fois"][0], "wb") as writer:
			writer.write("\n".join(foi))
	run_hypergeom(args["fois"][0], args["gfs"][0], args["bg_path"][0], args["run_files_dir"], "", False, "",
				  args['data_dir'], args["run_annotation"], run_randomization_test=False, organism=args['organism'],
				  stat_test=args['stat_test'])


if __name__ == "__main__":
	main()


