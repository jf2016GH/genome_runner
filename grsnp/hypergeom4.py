#!/usr/bin/env python2

import argparse
import collections
import math
import sys
import logging
from logging import FileHandler,StreamHandler
from bx.intervals.intersection import IntervalTree
from scipy.stats import hypergeom
import numpy
import scipy
import pdb
import os
import json
import rpy2.robjects as robjects
import gzip
import tarfile
import traceback
import StringIO
import dbcreator as bedfilecreator
import textwrap
import subprocess
import sys



# Logging configuration
logger = logging.getLogger()
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
# This line outputs logging info to the console


matrix_outpath = None
detailed_outpath = None
progress_outpath = None
output_dir = None
console_output = False 
print_progress = False
logger_path = "gr_log.txt"

def get_overlap_statistics(gf,fois):
    """Returns a dictionary with indicating how many hits exist for each foi against the gf
    gf: filepath for GF
    fois: list of FOI filepaths
    """
    results = []
    out = ""
    try:
        # Runs overlapStatistics with preprocessed background stats if they exist
        out = subprocess.Popen(["overlapStatistics"] + [gf] + fois,stdout=subprocess.PIPE)
        out.wait()
        tmp = out.stdout.read()
        if tmp[:6] == "ERROR:": 
            logger.error(tmp[7:])
            raise Exception(tmp)

        for x in tmp.split("\n")[1:]:
            if x != "":
                tmp = x.split("\t")
                foi_name,n,hit_count = os.path.split(tmp[0])[-1],tmp[2],tmp[3]
                results.append({"queryfile": foi_name,"queryregions": int(n),"intersectregions": int(hit_count),"indexregions": int(tmp[1])})
    except Exception, e:        
        logger.error(traceback.format_exc())
        return
    return results


def get_bgobs(bg,gf,bkg_overlap_path): 

    _write_progress("Getting overlap stats on background and {}".format(base_name(gf)))
    logger.info("Getting overlap stats on background and {}".format(base_name(gf)))
    # See if pre-calculated values exist
    if os.path.exists(bkg_overlap_path):
        data = open(bkg_overlap_path).read().split("\n")
        data = [x.split("\t") for x in data if x != ""]
        d_gf = [x[1] for x in data if x[0] == gf and x[1]  != ""]
        if len(d_gf) != 0:
            bg_obs = [x.split(":")[1] for x in d_gf[0].split(",") if x.split(":")[0] == bg]
            if len(bg_obs) != 0:
                logger.info("Pre-calculated values found for background and {} ".format(base_name(gf)))
                return bg_obs[0]


    # manually get overlap values
    result = get_overlap_statistics(gf,[bg])
    try:
        result = int(result[0]["intersectregions"])
    except Exception, e:
        result = None
        logger.error(traceback.format_exc())
    return result




def p_value(foi_obs,n_fois,bg_obs,n_bgs,foi_name,gf_name):    
    """Return the signed p-value of all FOIs against the GF.
    """
    global logger_path
    bg_obs,n_bgs = int(bg_obs),int(n_bgs)
    ctable = [[foi_obs, n_fois-foi_obs],
              [bg_obs-foi_obs,n_bgs-n_fois-(bg_obs-foi_obs)]]
              
    # Ensure there are no negative values in the ctable
    do_chi_square = True
    for i in ctable:
        for k in i:
            if k < 0:
                logger.warning("Cannot calculate p-value for {} and {}. Is the background too small?".format(gf_name,foi_name))
                return 1
            if k < 11:
                do_chi_square = False

    if bg_obs < foi_obs:
        odds_ratio, pval = "nan", 1
        logger.error("P-value cannot be calculated (pvalue = 1.0, odds_ratio = 'nan'). Number of SNPs overlapping with GF > number of background SNPs overlapping with GF.".format(foi_name,gf_name))
    else: 
        if do_chi_square:        
            chi_result = scipy.stats.chi2_contingency(ctable)
            odds_ratio = (ctable[0][0]*ctable[1][1])/(ctable[0][1]*ctable[1][0])
            pval = chi_result[1]
        else:    
            odds_ratio, pval = scipy.stats.fisher_exact(ctable)
    sign = -1 if (odds_ratio < 1) else 1

    # write out to the detailed results file
    write_output("\t".join(map(str, [foi_name.rpartition('/')[-1], foi_obs, n_fois, bg_obs, n_bgs, 
                "%.2f" % odds_ratio if type(odds_ratio) != type("") else odds_ratio, 
                "%.2e" % pval if type(pval) != type("") else pval,
                "Chi-squaredquared" if do_chi_square else "Fisher-Exact"])) + "\n",detailed_outpath)

    # write out to the enrichment result file
    er_result_path = os.path.join(output_dir,"enrichment")
    if not os.path.exists(er_result_path): os.mkdir(er_result_path)
    er_result_path = os.path.join(er_result_path,base_name(foi_name)+".txt")
    # writes the first line as the header line
    if not os.path.exists(er_result_path): write_output(foi_name+"\tP-value\tDirection\n",er_result_path)
    if sign == 1 or str(odds_ratio) == "inf":
        direction  = "overrepresented"    
    else: direction =  "underrepresented"
    if pval > 0.05:
        direction = "nonsignificant"
    write_output("\t".join([gf_name,"%.2e" % pval if type(pval) != type("") else pval,direction])+"\n",er_result_path) 
    if pval < 1E-307:
        # set to value obtained from sys.float_info.min_10_exp
        pval = 1E-306   
    return sign * pval




def cluster_matrix(input_path,output_path):
    '''
    Takes the matrix file outputted by genomerunner and clusters it in R
    '''        
    pdf_outpath = ".".join(output_path.split(".")[:-1] + ["pdf"])
    r_script = """library(gplots)
                    t5 = as.matrix(read.table("{}")) 
                    t5<-as.matrix(t5[apply(t5, 1, function(row) {{sum(abs(row) < 0.05) >= 1}}), ]) # Remove rows with all values below cutoff 2 (p-value 0.01)
                    if (nrow(t5) > 0 && ncol(t5) > 0) {{
                        # Log transform matrix and keep correct sign
                        for (i in 1:nrow(t5)) {{
                          for (j in 1:ncol(t5)) {{
                            if (t5[i, j] < 0) {{ t5[i, j] <- log10(abs(t5[i, j]))}} else {{ t5[i,j] <- -log10(t5[i,j])}}
                          }}
                        }}
                        if (dim(t5)[1] > 1 && dim(t5)[2] > 1) {{
                           pdf(file="{}")
                           h = heatmap.2(t5,  hclustfun=function(e) hclust(e,method="average"),margins=c(15,15),  distfun=function(x) dist(x,method="euclidean"), cexCol=1, cexRow=1)
                           dev.off()
                           write.table(t(h$carpet),"{}",sep="\t")
                        }} else {{
                            write.table(paste("ERROR: Cannot run clustering on",nrow(t5),"x",ncol(t5),"matrix. Should be at least 2 x 2. Analyze more sets of SNPs and select more genomic features"),"{}",sep="\t", row.names=F, col.names=F)
                        }}
                    }} else {{
                        write.table("ERROR: Nothing significant","{}",sep="\t",row.names=F,col.names=F)
                    }}""".format(input_path,pdf_outpath,output_path,output_path,output_path)
    robjects.r(r_script)
    return output_path    

def pearsons_cor_matrix(matrix_path,out_dir):
    global logger_path
    output_path = os.path.join(out_dir,"pcc_matrix.txt")
    pdf_outpath = ".".join(output_path.split(".")[:-1] + [".pdf"])
   
    #pcc = open(output_path).read() 
    #if "PCC can't be performed" not in pcc:
    r_script = """t5 = as.matrix(read.table(\""""+matrix_path+"""\")) 
        library(Hmisc)
        library(gplots)
        t5<-as.matrix(t5[,apply(t5,2,sd)!=0]) # Remove columns with SD = zeros
        if (dim(t5)[1] > 4 && dim(t5)[2] > 1) {
            p5<-rcorr(t5)
            pdf(file=\"""" + pdf_outpath +"""\")
            h<-heatmap.2(as.matrix(p5[[1]]),margins=c(15,15)) # [[1]] element contains actual PCCs, we cluster them
            dev.off()
            write.table(h$carpet,\"""" + output_path + """\",sep="\t") # Write clustering results
            write.table(p5[[3]][h$rowInd, h$colInd],\""""+output_path.split(".")[0]+"_pvalue.txt"+ """\",sep="\t") # [[3]] element contains p-values. We write them using clustering order
        } else {
            write.table(paste("ERROR: Cannot run correlation analysis on", dim(t5)[1], "x", dim(t5)[2], "matrix. Should be at least 5 x 2. Analyze more sets of SNPs and select more genomic features"),\"""" + output_path + """\",sep="\t", row.names=F, col.names=F) # Write clustering results            
        }"""
    robjects.r(r_script)
    return output_path


def get_annotation(foi,gfs):
    """
    fois: list of FOI filepath
    gfs: filepaths for GF
    """
    results = []
    out = subprocess.Popen(["annotationAnalysis"] + [foi] + gfs,stdout=subprocess.PIPE) # TODO enable ["--print-region-name"]
    out.wait()
    return out.stdout.read()

# Writes the output to the file specified.  Also prints to console if console_output is set to true
def write_output(content,outpath=None):
    if outpath:
        with open(outpath,"a") as out:
            out.write(content)
    if console_output:
        print >> sys.stderr, content

# Collect lines from a file into a list
def read_lines(path):
    elems = []
    with open(path) as h:
        for line in h:
            if line.strip():
                elems.append(line.strip())
    return elems

def base_name(k):
    return os.path.basename(k).split(".")[0]

def _write_progress(line):
    """Saves the current progress to the progress file
    """
    global progress_outpath
    if progress_outpath:
        global curprog, progmax
        progress = {"status": line, "curprog": curprog,"progmax": progmax}
        with open(progress_outpath,"wb") as progfile:
            progfile.write(json.dumps(progress))
    if print_progress:
        print line


def _write_head(content,outpath):
    f = front_appender(outpath)
    f.write(content)
    f.close()


def check_background_foi_overlap(bg,fois):
    """ Calculates the overlap of the FOIs with the background.
    Removes FOIs that are poorly formed with the background.
    """
    good_fois = []
    # Runs overlapStatistics on background and FOIs
    foi_bg_stats =  get_overlap_statistics(bg,fois)
    for f in foi_bg_stats:
        isgood = True
        foi_name,n_bgs,n_fois,foi_in = f["queryfile"],f["indexregions"],f["queryregions"],f["intersectregions"]
        if n_fois < 5:
            isgood = False
            logger.error("Number of SNPs in {} < 5. Removing it from analysis.".format(foi_name))
        elif n_bgs < n_fois:
            isgood = False
            logger.error("Number of SNPs in {} > than in background. Removing it from analysis.".format(foi_name))
        if isgood:
            # ensure that overlapStatistics output filename with extension for queryFile field
            good_fois.append([x for x in fois if os.path.basename(x) == f["queryfile"]][0])
        if foi_in < n_fois:
            logger.error("{} out of {} {} SNPs are not a part of the background. P-value are unreliable. Please, include all SNPs in the background and re-run analysis.".format(n_fois-foi_in,n_fois,foi_name))
    return [foi_bg_stats, good_fois]
                                                                                                                       


def get_description(gf,trackdb):
    desc = [x["longLabel"] for x in trackdb if x["tableName"] == gf]
    if len(desc) is not 0: return desc[0]
    else: return "No Description"




def _zip_run_files(fois,gfs,bg_path,outdir,id=""):
    '''
    File paths of FOIs and GFs as a list. Gathers all the files together in one zipped file
    '''    
    f = open(os.path.join(outdir,"gr_log.txt"))
    f_log = f.read()
    f.close()
    path_settings,f_sett =os.path.join(outdir,".settings"),""
    if os.path.exists(path_settings):
        f = open(path_settings)
        f_sett = f.read() + "\n###LOG###\n"
        f.close()
    new_log_path = os.path.join(outdir,"gr_log.txt")
    new_log = open(new_log_path,'wb')
    new_log.write(f_sett+f_log)
    new_log.close()
    tar_path = os.path.join(outdir,'GR_{}.tar'.format(id))
    tar = tarfile.TarFile(tar_path,"a")    
    output_files =  [os.path.join(outdir,x) for x in os.listdir(outdir) if x.endswith(".txt") or x.endswith(".pdf")]
    fls = output_files
    for f in fls:
        tar.add(f,os.path.basename(f))
    tar.close()
    tar_file = open(tar_path,'rb')
    with gzip.open(tar_path+".gz","wb") as gz:
        gz.writelines(tar_file)
    tar_file.close()
    if os.path.exists(tar_path): os.remove(tar_path)

def run_hypergeom(fois, gfs, bg_path,outdir,job_name="",zip_run_files=False,bkg_overlaps_path="",gr_data_dir = "" ,run_annotation=False):
    global formatter
    global detailed_outpath,matrix_outpath, progress_outpath, curprog, progmax,output_dir
    if not os.path.exists(os.path.normpath(outdir)): os.mkdir(os.path.normpath(outdir))
    sett_path = os.path.join(outdir,".settings")
    fh = logging.FileHandler(os.path.join(outdir,'gr_log.txt'))
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    logger.setLevel(logging.INFO)
    output_dir = outdir
    curprog,progmax = 0,1
    organism = ""
    try:
        trackdb = []
        # loads the trackDb file data so that descriptions for the GF can be outputted
        if os.path.exists(sett_path): 
            with open(sett_path) as re:
                organism = [x.split("\t")[1] for x in re.read().split("\n") if x.split("\t")[0] == "Organism:"][0]  # read in th organism data from .settings file generated by the server
                trackdb_path = os.path.join(gr_data_dir,organism,"trackDb")
                if os.path.exists(trackdb_path+".txt.gz") and os.path.exists(trackdb_path + ".sql"):
                    trackdb = bedfilecreator.load_tabledata_dumpfiles(trackdb_path)
        # set output settings
        detailed_outpath =  os.path.join(outdir, "detailed.txt") 
        matrix_outpath = os.path.join(outdir,"matrix.txt")
        progress_outpath = os.path.join(outdir,".prog")
        f = open(matrix_outpath,'wb') 
        f.close()
        f = open(detailed_outpath,'wb')
        f.close()
        hdlr = logging.FileHandler(logger_path)
        formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
        hdlr.setFormatter(formatter)
        logger.addHandler(hdlr)
        logger.propagate = False

        # Read in the paths
        fois = [line for line in read_lines(fois) if not line.endswith(".tbi")]
        gfs = [line for line in read_lines(gfs) if not line.endswith(".tbi")]
        if bg_path.endswith(".tbi"):
            logger.error("Background has invalid extension (.tbi). Terminating run.")
            _write_progress("ERROR: Background has invalid extension (.tbi). Terminating run.")
            return

        foi_bg,good_fois = check_background_foi_overlap(bg_path,fois)   # Validate FOIs against background. Also get the size of the background (n_bgs)
        write_output("\t".join(map(base_name,good_fois))+"\n", matrix_outpath)
        write_output("\t".join(['foi_name', 'foi_obs', 'n_fois', 'bg_obs', 'n_bgs', 'odds_ratio', 'p_val','test_type']) + "\n",detailed_outpath)
        curprog,progmax = 0,len(gfs)
        _write_progress("Performing calculations on the background.")
        for gf in gfs: 
            current_gf = base_name(gf)      
            _write_progress("Performing Hypergeometric analysis for {}".format(base_name(gf)))
            write_output("###"+base_name(gf)+"\t"+get_description(base_name(gf),trackdb)+"###"+"\n",detailed_outpath)
            res = get_overlap_statistics(gf,good_fois) 

            # calculate bg_obs
            bg_obs = get_bgobs(bg_path,gf,bkg_overlaps_path)
            if bg_obs == None: 
                logger.error("Skipping {}".format(gf))
                continue

            n_bgs = foi_bg[0]["indexregions"]
            # run the enrichment analysis and output the matrix line for the current gf
            write_output("\t".join([base_name(gf)] + [str(p_value(res[i]["intersectregions"],res[i]["queryregions"],bg_obs,n_bgs ,os.path.basename(good_fois[i]),os.path.basename(gf))) for i in range(len(good_fois))])+"\n",matrix_outpath)
            curprog += 1
        if len(gfs) > 1 and len(good_fois) > 1:
            clust_path =  cluster_matrix(matrix_outpath,os.path.join(outdir,"clustered.txt"))
            if len(gfs) > 4:               
                pearsons_cor_matrix(clust_path,outdir)
            else:
                with open(os.path.join(os.path.join(outdir,"pcc_matrix.txt")),"wb") as wb:
                    wb.write("ERROR:PCC matrix requires at least a 5 X 2 matrix.")
        else:
            with open(os.path.join(outdir,"clustered.txt"),"wb") as wb:
                wb.write("ERROR:Clustered matrix requires at least a 2 X 2 matrix.")
        if run_annotation:
            annot_outdir = os.path.join(outdir,"annotations")
            if not os.path.exists(annot_outdir): os.mkdir(annot_outdir)
            curprog,progmax = 0,len(fois)
            for f in fois:                
                _write_progress("Running Annotation Analysis for {}.".format(base_name(f)))
                with open(os.path.join(annot_outdir,base_name(f) + ".txt"),"wb") as wr:
                    anot = get_annotation(f,gfs).split("\n")
                    wr.write("\t".join(base_name(x) for x in anot[0].split("\t")) + "\tTotal")
                    ind = 1
                    for a in anot[1:]:
                        if a != "":
                            cur_row = a.split("\t")
                            wr.write("\n" + str(ind) + "|"+"\t".join(cur_row + [str(sum([int(x) for x in cur_row[1:] if x != ""]))]))
                            ind += 1
                curprog += 1

        if zip_run_files:
            _write_progress("Preparing run files for download")
            _zip_run_files(fois,gfs,bg_path,outdir,job_name)
        _write_progress("Analysis Completed")       
    except Exception, e: 
        logger.error( traceback.print_exc())
        write_output(traceback.format_exc(),logger_path)
        _write_progress("Run crashed. See end of log for details.")

if __name__ == "__main__":
    global print_progress
    print_progress = True
    parser = argparse.ArgumentParser(description="Analysis of several FOI files against several GFs using Fisher's exact test or Chi-square. Best used for SNP set analysis, using whole SNP database as a spot background.")
    parser.add_argument("fois", nargs=1, help="Text file with FOI file names (SNPs only).") 
    parser.add_argument("gfs" ,nargs=1, help="Text file with GF file names, gzipped.") 
    parser.add_argument("bg_path", nargs=1, help="Path to spot background file (SNPs only).")
    parser.add_argument("--run_annotation" , "-a", help="Run annotation analysis", action="store_true" )
    parser.add_argument("--output_dir","-d", help="Directory to output the result to", default="")
    args = vars(parser.parse_args())
    run_hypergeom(args["fois"][0],args["gfs"][0],args["bg_path"][0],args["output_dir"],"",False,"","",args["run_annotation"])



class front_appender:
    '''
    Appends content to start of file.
    '''
    def __init__(self, fname, mode='a'):
        self.__write_queue = []
        self.__old_content = ""
        if mode == 'a': 
            self.__old_content = open(fname).read()
        self.__f = open(fname, 'w')


    def write(self, s):
        self.__write_queue.append(s)

    def close(self):
        self.__f.writelines(self.__write_queue + [self.__old_content])
        self.__f.close()

