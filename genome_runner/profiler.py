import os
import cProfile
import query
import argparse


def runprofile(fois,gfs):
	for f in fois:
		query.run_enrichments(1000000,f,gfs,3,None,None,None,'hg19')

def get_profile_outpath():
	n = 0
	outputpath = os.path.join("profile","{}.profile".format(n))
	while os.path.exists(outputpath): 
		n+=1
		outputpath = os.path.join("profile",str(n) + ".profile")
	return outputpath


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="""Profile the GenomeRunner Enrichment Analysis, WARNING: do not run in production environment, 
			will cause enrichment results with the same name to be overwritten""")
	parser.add_argument('--output','-o', help='The folder to output the profile data to')
	parser.add_argument('--inputfoi','-i', help='The folder that contains the foi files')
	parser.add_argument('--inputgf','-g',help='The folder that contains the foi files')
	args = vars(parser.parse_args())
	foidir = os.path.join("profile","foi")
	gfdir = os.path.join("profile","gf")
	gf = [os.path.join(gfdir,each) for each in os.listdir(gfdir) if each.endswith(('.bed','.gz') )]
	print "The following files will be considered as genomic feature files.\n {}".format(gf)
	foi = [os.path.join(foidir,each) for each in os.listdir(foidir) if each.endswith(('.bed','.gz'))]
	print "The following files will be consdered as features of interest files. \n {}".format(foi)
	profileoutpath = get_profile_outpath()
	print profileoutpath
	cProfile.run("runprofile(gf,foi)",profileoutpath)
	print "profiling done! Outputed to {}".format(profileoutpath)


