#!/usr/bin/env python2
"""
Analysis of several FOI files against several GFs using Fisher's exact test. Best used for SNP set analysis, using whole SNP database as a spot background.
sys.argv[1] - text file with FOI file names. Include full, or relative path, if needed.
sys.argv[2] - text file with GF file names. Include full, or relative path, if needed.
sys.argv[3] - spot background file
"""

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
import bedfilecreator
import textwrap
import subprocess
import sys



# Logging configuration
logger = logging.getLogger()
logger = logging.getLogger('genomerunner.query')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
# This line outputs logging info to the console
logger.setLevel(logging.INFO)


matrix_outpath = None
detailed_outpath = None
progress_outpath = None
console_output = False 
logger_path = "log.txt"

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

    if os.path.exists(bkg_overlap_path):
        _write_progress("Getting overlap stats on background and {}".format(gf))
        write_output("Getting overlap stats on background and {}".format(gf),logger_path)
        data = open(bkg_overlap_path).read().split("\n")
        data = [x.split("\t") for x in data if x != ""]
        d_gf = [x[1] for x in data if x[0] == gf and x[1]  != ""]
        if len(d_gf) != 0:
            bg_obs = [x.split(":")[1] for x in d_gf[0].split(",") if x.split(":")[0] == bg]
            if len(bg_obs) != 0:
                return bg_obs[0]
    logger.info("Manually calculating background stats")
    result = get_overlap_statistics(bg,[gf])
    return int(result[0]["intersectregions"])




def p_value(foi_obs,n_fois,bg_obs,n_bgs,foi_name,gf_name):    
    """Return the signed log10 p-value of all FOIs against the GF.
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
                return math.log10(1)
            if k < 11:
                do_chi_square = False



    if n_fois == foi_obs or n_bgs == bg_obs: 
        odds_ratio, pval = "nan", 1
        logger.warning("P-value cannot be calculated (set to 1.0). Number of {} equal to # SNPs overlapping with {} (GF) or number of SNPs in background equal to number overlapping with GF".format(foi_name,gf_name))
    elif n_fois < 5: 
        odds_ratio, pval = "nan", 1
        logger.warning("P-value cannot be calculated for {}, must have > 5 (p-value set to 1.0)".format(foi_name))
    else: 
        if do_chi_square:        
            logger.info("Using the Chi-squared test for {} and {}. Ctable values all > 10: {}".format(gf_name,foi_name,ctable))
            chi_result = scipy.stats.chi2_contingency(ctable)
            obs = ctable[0][0]
            exp = chi_result[3][0][0]
            odds_ratio = obs/exp
            pval = chi_result[1]
        else:    
            logger.info("Using the Fisher's exact test test for {} and {}. Ctable values not all > 10: {}".format(gf_name,foi_name,ctable))
            odds_ratio, pval = scipy.stats.fisher_exact(ctable)
    sign = 1 if (odds_ratio < 1) else -1
    write_output("\t".join(map(str, [foi_name.rpartition('/')[-1], foi_obs, n_fois, bg_obs, n_bgs, 
                "%.2f" % odds_ratio if type(odds_ratio) != type("") else odds_ratio, 
                "%.2f" % pval if type(pval) != type("") else pval])) + "\n",detailed_outpath)  
    
    return sign * math.log10(pval)   




def cluster_matrix(input_path,output_path):
    '''
    Takes the matrix file outputted by genomerunner and clusters it in R
    '''        
    r_script = """library(gplots)
                    t5 = read.table("{}")                    
                    pt5 <- as.matrix(t5)
                    row_names <- rownames(pt5)
                    col_names <- colnames(pt5)
                    pt5 <- apply(pt5, 2, function(x) as.numeric(x))
                    row.names(pt5) <- row_names
                    colnames(pt5) <- col_names
                    if (!all(pt5[1,1] == pt5))
                    {{
                       h = heatmap.2(pt5,  hclustfun=function(m) hclust(m,method="average"),  distfun=function(x) dist(x,method="euclidean"), cexCol=1, cexRow=1)
                       write.table(t(h$carpet),"{}",sep="\t")
                    }} else {{
                        write.table(t5,"{}",sep="\t")
                    }}""".format(input_path,output_path,output_path)
    robjects.r(r_script)
    return output_path    

def pearsons_cor_matrix(matrix_path,out_dir):
    global logger_path
    output_path = os.path.join(out_dir,"pcc_matrix.txt")
    outputpdf_path = ".".join(output_path.split(".")[:-1] + [".pdf"])
    with open(matrix_path) as f:
        # check if each row has unique values
        rows = [x.rstrip() for x in f if x != ""]
        unique = True
        for k in rows[1:]:
            if len(k) is 0: continue
            lst = k.split("\t")[1:]
            if lst[1:] == lst[:-1]: # check if all values in the list are the same
                unique = False   

        # check if each rolumn has unique values
        len_row = len(rows[0].split("\t"))
        for k in range(1,len_row):
            lst = [rows[x].split("\t")[k] for x in range(1,len(rows))]  
            if lst[1:] == lst[:-1]: 
                unique = False   

    if unique == False:
        with open(output_path, "wb") as f:  
            f.write("ERROR:PCC cannot be performed.  Each row/column of matrix must have varying values.")
        logger.warning("PCC cannot be performed.  Each row/column of matrix must have varying values.")
        return output_path

    ### this calculates the PCC matrix
    r_script = """library(gplots)
                    t5 = read.table("{}")                    
                    pt5 <- cor(as.matrix(t5)) # Correlation matrix
                    row_names <- rownames(pt5)
                    col_names <- colnames(pt5)
                    pt5 <- apply(pt5, 2, function(x) as.numeric(x))
                    row.names(pt5) <- row_names
                    colnames(pt5) <- col_names 
                    h = heatmap.2(pt5,  hclustfun=function(m) hclust(m,method="average"),  distfun=function(x) dist(x,method="euclidean"), cexCol=1, cexRow=1)
                    write.table(t(h$carpet),"{}",sep="\t")
                    """.format(matrix_path,outputpdf_path,output_path)
    #robjects.r(r_script)

    ## calculates the PCC p-values 

    r_script = """t5 = read.table(\""""+output_path+"""\")                   
                
                  pn <- function(X){crossprod(!is.na(X))}

                  col <- function (x, as.factor = FALSE) 
                    {
                      if (as.factor) {
                        labs <- colnames(x, do.NULL = FALSE, prefix = "")
                        res <- factor(.Internal(col(dim(x))), labels = labs)
                        dim(res) <- dim(x)
                        res
                      }
                      else .Internal(col(dim(x)))
                    }
                cor.prob <- function(X){
                  pair.SampSize <- pn(X)
                  above1 <- row(pair.SampSize) < col(pair.SampSize)
                  pair.df <- pair.SampSize[above1] - 2
                  R <- cor(X, use="pair")
                  above2 <- row(R) < col(R)
                  r2 <- R[above2]^2
                  Fstat <- (r2 * pair.df)/(1 - r2)
                  R[above2] <- 1 - pf(Fstat, 1, pair.df)
                  R
                }

                pt5.prob <- cor.prob(as.matrix(t5))
                row_names <- rownames(pt5)
                col_names <- colnames(pt5)
                row.names(pt5.prob) <- row_names
                colnames(pt5.prob) <- col_names 
                write.table(pt5.prob,\""""+output_path.split(".")[0]+"_pvalue.txt"+"""\",sep="\t")"""
    #robjects.r(r_script)
    r_script = """t5 = read.table(\""""+matrix_path+"""\") 
        library("Hmisc")
        p5<-rcorr(as.matrix(t5))
        p5[[1]]
        p5[[3]]
        h<-heatmap.2(as.matrix(p5[[1]]))
        write.table(h$carpet,\"""" + output_path + """\",sep="\t")
        write.table(p5[[3]][h$rowInd, h$colInd],\""""+output_path.split(".")[0]+"_pvalue.txt"+ """\",sep="\t")""" 
    robjects.r(r_script)
    return output_path   


def get_annotation(foi,gfs):
    """
    fois: list of FOI filepath
    gfs: filepaths for GF
    """
    results = []
    out = subprocess.Popen(["annotationAnalysis"] + ["--print-region-name"] + [foi] + gfs,stdout=subprocess.PIPE)
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

def base_name(path):
    return ".".join(os.path.basename(path).split(".")[:-1])

def _write_progress(line):
    """Saves the current progress to the progress file
    """
    global progress_outpath
    if progress_outpath:
        global curprog, progmax
        progress = {"status": line, "curprog": curprog,"progmax": progmax}
        with open(progress_outpath,"wb") as progfile:
            progfile.write(json.dumps(progress))


def _write_head(content,outpath):
    f = front_appender(outpath)
    f.write(content)
    f.close()


def check_background_foi_overlap(bg,fois):
    """ Calculates the overlap of the FOIs with the background
    """
    foi_bg_stats =  get_overlap_statistics(bg,fois)
    logger.info("###Background and SNPs stats###\n")


    return foi_bg_stats


def run_hypergeom(fois, gfs, bg_path,outdir,job_name="",zip_run_files=False,run_annotation=False):
    sett_path = os.path.join(outdir,".settings")
    logger_path = os.path.join(outdir,'log.txt')
    global detailed_outpath,matrix_outpath, progress_outpath, curprog, progmax
    curprog,progmax = 0,1
    try:
        trackdb = []
        if os.path.exists(sett_path):        
            with open(sett_path) as re:
                organism = [x.split("\t")[1] for x in re.read().split("\n") if x.split("\t")[0] == "Organism:"][0]
                trackdb_path = os.path.join("data",organism,"trackDb")
                if os.path.exists(trackdb_path):
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

        fois = read_lines(fois)
        gfs = read_lines(gfs)

        write_output("\t".join(map(base_name,fois))+"\n", matrix_outpath)
        write_output("\t".join(['foi_name', 'foi_obs', 'n_fois', 'bg_obs', 'n_bgs', 'odds_ratio', 'p_val']) + "\n",detailed_outpath)
        curprog,progmax = 0,len(gfs)
        _write_progress("Performing calculations on the background.")
        foi_bg = check_background_foi_overlap(bg_path,fois)
        for gf in gfs: 
            current_gf = base_name(gf)      
            _write_progress("Performing Hypergeometric analysis for {}".format(base_name(gf)))   
            write_output("###"+base_name(gf)+"\t"+get_description(base_name(gf),trackdb)+"###"+"\n",detailed_outpath)

            res = get_overlap_statistics(gf,fois)  
            bg_obs = get_bgobs(bg_path,gf,os.path.join("data",organism,"bkg_overlaps.gr"))
            # run the enrichment analysis and output the matrix line for the current gf
            write_output("\t".join([base_name(gf)] + [str(p_value(res[i]["intersectregions"],res[i]["queryregions"],bg_obs,foi_bg[0]["indexregions"] ,os.path.basename(fois[i]),os.path.basename(gf))) for i in range(len(fois))])+"\n",matrix_outpath)
            curprog += 1
        if len(gfs) > 1 and len(fois) > 1:
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
                            wr.write("\n" + str(ind) + "|"+"\t".join(cur_row + [str(sum([int(x) for x in cur_row[1:]]))]))
                            ind += 1
                curprog += 1

            
        _write_progress("Preparing run files for download")
        _zip_run_files(fois,gfs,bg_path,outdir,job_name)
        _write_progress("Analysis Completed")       
    except Exception, e: 
        logger.error( traceback.print_exc())
        write_output(traceback.format_exc(),logger_path)
        _write_progress("Run crashed. See end of log for details.")

def get_description(gf,trackdb):
    desc = [x["longLabel"] for x in trackdb if x["tableName"] == gf]
    if len(desc) is not 0: return desc[0]
    else: return "No Description"




def _zip_run_files(fois,gfs,bg_path,outdir,job_name=""):
    '''
    File paths of FOIs and GFs as a list. Gathers all the files together in one zipped file
    '''    
    f = open(os.path.join(outdir,"log.txt"))
    f_log = f.read()
    f.close()
    path_settings,f_sett =os.path.join(outdir,".settings"),""
    if os.path.exists(path_settings):
        f = open(path_settings)
        f_sett = f.read() + "\n###LOG###\n"
        f.close()
    new_log_path = os.path.join(outdir,"log.txt")
    new_log = open(new_log_path,'wb')
    new_log.write(f_sett+f_log)
    new_log.close()
    tar_path = os.path.join(outdir,'GR_Runfiles_{}.tar'.format(job_name))
    tar = tarfile.TarFile(tar_path,"a")    
    output_files =  [os.path.join(outdir,x) for x in os.listdir(outdir) if x.endswith(".txt")]
    fls = output_files
    for f in fls:
        tar.add(f,os.path.basename(f))
    tar.close()
    tar_file = open(tar_path,'rb')
    with gzip.open(tar_path+".gz","wb") as gz:
        gz.writelines(tar_file)
    tar_file.close()
    if os.path.exists(tar_path): os.remove(tar_path)


if __name__ == "__main__":
    console_output,logger_path = True,os.path.join(outdir,'.log')
    parser = argparse.ArgumentParser(description="Create a matrix of hypergeometric p-values for genomic intersections.")
    parser.add_argument("--fois", "-f", nargs=1, help="Text file with FOI file names (SNPs only).") 
    parser.add_argument("--gfs" , "-g",nargs=1, help="Text file with GF file names, gzipped.") 
    parser.add_argument("--bg_path" "-b", nargs=1, help="Path to spot background file (SNPs only).")
    parser.add_argument("--run_annotation" , "-a", help="Run annotation analysis", action="store_true" )

    args = parser.parse_args()
    fois = read_lines(args.fois[0])
    gfs = read_lines(args.gfs[0])

    bg = read_intervals(args.bg_path[0])
    bg_iset = IntervalSet(bg)    
    f = open(matrix_outpath,'wb')
    f.close()
    f = open(detailed_outpath,'wb')
    f.close()

    _write_head("\n\n#Detailed log report#\n",logger_path)
    foi_sets = dict((path,read_intervals(path, snp_only=True, background=bg_iset,d_path=logger_path)) for path in fois)
    _write_head("#Grooming Summary#",logger_path)

    write_output("\t".join(fois)+"\n", matrix_outpath)
    write_output("\t".join(['foi_name', 'foi_obs', 'n_fois', 'bg_obs', 'n_bgs', 'odds_ratio', 'p_val'])+"\n",detailed_outpath)
    global current_gf
    for gf in gfs:
        current_gf = base_name(gf)
        gf_iset = IntervalSet(read_intervals(gf))
        write_output(gf+"\n",detailed_outpath) 
        write_output("\t".join([gf] + [str(p_value(gf_iset, foi_sets[foi], bg, foi,True)) for foi in fois])+"\n",matrix_outpath)
    if len(gfs) > 1 and len(fois) > 1:
        cluster_matrix(matrix_outpath,os.path.join(outpath,"matrix_clustered.gr"))
        pearsons_cor_matrix(os.path.join(outdir,"matrix_clustered.gr"),outdir)
    else:
        with open(os.path.join(outpath,"matrix_clustered.gr"),"wb") as wb:
            wb.write("ERROR:Clustered matrix requires at least a 2 X 2 matrix.")


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

