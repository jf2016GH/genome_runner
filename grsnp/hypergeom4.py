#!/usr/bin/env python2
from __future__ import division
import argparse
import collections
import math
import sys
import logging
from logging import FileHandler,StreamHandler
#from bx.intervals.intersection import IntervalTree
from scipy.stats import hypergeom
import numpy as np
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
import commands
import mako
import simplejson
import zipfile
import inspect
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import grsnp.dbcreator_util as grsnp_util

# Logging configuration
logger = logging.getLogger()
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
logger.propagate = 0

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
        out = subprocess.Popen(["overlapStatistics"] + [gf] + fois,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out.wait()
        tmp = out.stdout.read()
	tmp_er = out.stderr.read()
	if tmp_er != "": logger.error(tmp_er)
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

def p_value(foi_obs,n_fois,bg_obs,n_bgs,foi_path,gf_path,background_path,run_randomization_test=False):    
    """Return the signed p-value of all FOIs against the GF.
    """
    global logger_path
    foi_name = base_name(foi_path)
    gf_name = base_name(gf_path)
    sign,pval,odds_ratio,did_chi_square = calculate_p_value(foi_obs,n_fois,bg_obs,n_bgs,foi_name,gf_path)


    # write out to the enrichment result file
    er_result_path = os.path.join(output_dir,"enrichment")
    if not os.path.exists(er_result_path): os.mkdir(er_result_path)
    er_result_path = os.path.join(er_result_path,base_name(foi_name)+".txt")
    # writes the first line as the header line
    if not os.path.exists(er_result_path): write_output("GF Name\tP-value\tDirection\n",er_result_path)
    if sign == 1 or str(odds_ratio) == "inf":
        direction  = "overrepresented" 
    else: direction =  "underrepresented"

    # calculate the p_rand
    prnd = 1 # default prnd for non-significant results
    if pval > 0.05:
        direction = "nonsignificant"
    else:
        if run_randomization_test: 
            _write_progress("Running randomization test on {}".format(foi_name))
            prnd = p_rand(foi_path,n_fois,background_path,bg_obs,n_bgs,gf_path)  
    
    pval_unmod = pval
    pval = np.power(10,-(np.log10(prnd)- np.log10(pval))) # adjust p_value using randomization test
    # write out to the detailed results file
    strpval,strprnd = "","" 
    if run_randomization_test:
        strpval = "%.2e" % pval if type(pval) != type("") else pval 
        strprnd = "%.2e" % prnd if type(prnd) != type("") else prnd 
    write_output("\t".join(map(str, [foi_name.rpartition('/')[-1], foi_obs, n_fois, bg_obs, n_bgs, 
                "%.2f" % odds_ratio if type(odds_ratio) != type("") else odds_ratio, 
                "%.2e" % pval_unmod if type(pval_unmod) != type("") else pval_unmod,
                "Chi-squared" if did_chi_square else "Fisher-Exact",
                strprnd,strpval])) + "\n",detailed_outpath)

    write_output("\t".join([gf_name,"%.2e" % pval if type(pval) != type("") else pval,direction])+"\n",er_result_path) 
    if pval < 1E-307:
        # set to value obtained from sys.float_info.min_10_exp
        pval = 1E-306   
    return sign * pval

def p_rand(foi_path,n_fois,background_path,bg_obs,n_bgs,gf_path):
    ''' Calculated by generating 'num' random feature files and running them against gf_path.
    Calculates the mean of the p_values for the overrepresented and underrepresented random features separately.
    '''
    num = 10
    rnds_paths = generate_randomsnps(foi_path,background_path,n_fois,gf_path,num)
    rnd_stats = get_overlap_statistics(gf_path,rnds_paths)
    p_rand = [1]
    for r in rnd_stats:
        sign,pval,odds_ratio,chi = calculate_p_value(r["intersectregions"],r["queryregions"],bg_obs,n_bgs,base_name(foi_path),gf_path)
        p_rand.append(pval)
    return np.min(p_rand)

def calculate_p_value(foi_obs,n_fois,bg_obs,n_bgs,foi_name,gf_path):
    bg_obs,n_bgs = int(bg_obs),int(n_bgs)
    ctable = [[foi_obs, n_fois-foi_obs],
              [bg_obs-foi_obs,n_bgs-n_fois-(bg_obs-foi_obs)]]
              
    # Ensure there are no negative values in the ctable
    do_chi_square = True
    for i in ctable:
        for k in i:
            if k < 0:
                logger.warning("Cannot calculate p-value for {} and {}. Is the background too small? foi_obs {}, n_fois {}, bg_obs {}, n_bgs {}".format(base_name(gf_path),foi_name,foi_obs,n_fois,bg_obs,n_bgs))
                return [1,1,1,False]
            if k < 11:
                do_chi_square = False

    if bg_obs < foi_obs:
        odds_ratio, pval = "nan", 1
        logger.error("P-value cannot be calculated (pvalue = 1.0, odds_ratio = 'nan'). Number of SNPs overlapping with GF > number of background SNPs overlapping with GF. foi_obs {}, n_fois {}, bg_obs {}, n_bgs {}".format(foi_name,gf_name,foi_name,foi_obs,n_fois,bg_obs,n_bgs))
    else: 
        if do_chi_square:        
            chi_result = scipy.stats.chi2_contingency(ctable)
            odds_ratio = (ctable[0][0]*ctable[1][1])/(ctable[0][1]*ctable[1][0])
            pval = chi_result[1]
        else:    
            odds_ratio, pval = scipy.stats.fisher_exact(ctable)
    sign = -1 if (odds_ratio < 1) else 1
    return [sign,pval,odds_ratio,do_chi_square]






def generate_randomsnps(foi_path,background,n_fois,gf_path,num):
    paths = []
    out_dir = os.path.join(os.path.dirname(foi_path),"random")
    if not os.path.exists(out_dir): os.mkdir(out_dir)
    for n in range(num):        
        rnd_snp_path = os.path.join(out_dir,"random{}_".format(n)+base_name(foi_path)+".bed")
        # generate random snps from background
        if background.endswith('.gz'):
            command = "zcat {} | shuf -n {} | cut -f 1-3 > {}".format(background,str(n_fois),rnd_snp_path)
            out = commands.getstatusoutput(command)
        else:        
            command = "shuf -n {} {}  | cut -f 1-3 > {}".format(str(n_fois),str(background),rnd_snp_path)
            out = commands.getstatusoutput(command)
        paths.append(rnd_snp_path)
    return paths

def adjust_pvalue(input_path,method):
    ''' Uses R to adjust the pvalues
    '''
    try:
        if method == "None": return
        saved_stdout, saved_stderr = sys.stdout, sys.stderr
        sys.stdout = sys.stderr = open(os.devnull, "w")
        r_script = """
            t5 = as.matrix(read.table(\""""+input_path+"""\",sep="\t", header=T, row.names=1))

            rptemp <- t5

            rp1 <- apply(abs(t5), 2, function(x) {p.adjust(x, \""""+method+"""\")})

            for (i in 1:nrow(t5)){
              for (j in (1:ncol(t5))) {
                if (rptemp[i,j] < 0) { rp1[i,j] <- -1*rp1[i,j]}
              }
            }
            write.table(rp1,\""""+input_path+"""\",sep="\t")
        """
        robjects.r(r_script)
        sys.stdout, sys.stderr = saved_stdout, saved_stderr
        
    except Exception, e:
        logger.error("R CRASHED")
        logger.error(traceback.print_exc())
        logger.error(str(e))

def adjust_detailed_pvalue(input_path,method):   
    try:
        if method == "None": return
        saved_stdout, saved_stderr = sys.stdout, sys.stderr
        sys.stdout = sys.stderr = open(os.devnull, "w")
        r_script = """
            t5 <- as.data.frame(read.table(\""""+input_path+"""\", sep="\t", header=TRUE))

            rp1 <- p.adjust(abs(t5$P.value))
            for (i in 1:length(t5$P.value)) {
              if(t5$P.value[i] < 0) {
                rp1[i] <- -1*rp1[i]
              } 
            }

            t5 <- cbind(t5, P.adj=rp1)

            write.table(t5, \""""+input_path+"""\", sep="\t", row.names=F,quote=F)

        """
        robjects.r(r_script)
        sys.stdout, sys.stderr = saved_stdout, saved_stderr
        
    except Exception, e:
        logger.error("R CRASHED")
        logger.error(traceback.print_exc())        
        logger.error(str(e))

def cluster_matrix(input_path,output_path):
    '''
    Takes the matrix file outputted by genomerunner and clusters it in R
    '''        
    pdf_outpath = ".".join(output_path.split(".")[:-1] + ["pdf"])
    saved_stdout, saved_stderr = sys.stdout, sys.stderr
    sys.stdout = sys.stderr = open(os.devnull, "w")
    r_script = """t5 = as.matrix(read.table("{}"))                     
                    t5<-as.matrix(t5[apply(t5, 1, function(row) {{sum(abs(row) < 0.0001) >= 1}}), ]) # Remove rows with all values below cutoff 2 (p-value 0.01)
                    if (nrow(t5) > 0 && ncol(t5) > 0) {{
                        # Log transform matrix and keep correct sign
                        for (i in 1:nrow(t5)) {{
                          for (j in 1:ncol(t5)) {{
                            if (t5[i, j] < 0) {{ t5[i, j] <- log10(abs(t5[i, j]))}} else {{ t5[i,j] <- -log10(t5[i,j])}}
                          }}
                        }}
                        if (dim(t5)[1] > 1 && dim(t5)[2] > 1) {{
                            library(gplots)
                            library(RColorBrewer)
                            color<-colorRampPalette(c("blue","yellow"))
                            pdf(file="{}")
                            h = heatmap.2(t5, margins=c(20,20), col=color, trace="none", density.info="none", cexRow=1/log10(nrow(t5)),cexCol=1/log10(nrow(t5)))
                            dev.off()
                            write.table(t(h$carpet),"{}",sep="\t")
                        }} else {{
                            write.table(paste("ERROR: Cannot run clustering on",nrow(t5),"x",ncol(t5),"matrix. Should be at least 2 x 2. Analyze more sets of SNPs and select more genomic features"),"{}",sep="\t", row.names=F, col.names=F)
                        }}
                    }} else {{
                        write.table("ERROR: Nothing significant","{}",sep="\t",row.names=F,col.names=F)
                    }}""".format(input_path,pdf_outpath,output_path,output_path,output_path)
    robjects.r(r_script)
    sys.stdout, sys.stderr = saved_stdout, saved_stderr
    # write matrices in json format     
    json_mat = []       
    mat = open(output_path).read()      
    json_mat.append({"log": True,"neg": "Underrepresented", "pos": "Overrepresented","name": "P-value","alpha": 0.05,"matrix": mat})        
    with open(".".join(output_path.split(".")[:-1]+ ["json"]),'wb') as f:       
        simplejson.dump(json_mat,f)
    return output_path    

def pearsons_cor_matrix(matrix_path,out_dir):
    global logger_path
    output_path = os.path.join(out_dir,"pcc_matrix.txt")
    pval_output_path = os.path.join(os.path.split(output_path)[0], base_name(output_path)+"_pvalue.txt")
    pdf_outpath = ".".join(output_path.split(".")[:-1] + [".pdf"])
   
    saved_stdout, saved_stderr = sys.stdout, sys.stderr
    sys.stdout = sys.stderr = open(os.devnull, "w")
    #pcc = open(output_path).read() 
    #if "PCC can't be performed" not in pcc:
    r_script = """t5 = as.matrix(read.table(\""""+matrix_path+"""\")) 
        t5<-as.matrix(t5[,apply(t5,2,sd)!=0]) # Remove columns with SD = zeros
            if (dim(t5)[1] > 4 && dim(t5)[2] > 1) {
            library(Hmisc)
            library(gplots)
            library(RColorBrewer)
            color<-colorRampPalette(c("blue","yellow"))
            p5<-rcorr(t5)
            pdf(file=\"""" + pdf_outpath +"""\")
            h<-heatmap.2(as.matrix(p5[[1]]),margins=c(25,25), col=color, trace="none", density.info="none", cexRow=1/log10(nrow(p5[[1]])),cexCol=1/log10(nrow(p5[[1]]))) # [[1]] element contains actual PCCs, we cluster them
            dev.off()
            write.table(h$carpet,\"""" + output_path + """\",sep="\t") # Write clustering results
            write.table(p5[[3]][h$rowInd, h$colInd],\""""+pval_output_path+ """\",sep="\t") # [[3]] element contains p-values. We write them using clustering order
        } else {
            write.table(paste("ERROR: Cannot run correlation analysis on", dim(t5)[1], "x", dim(t5)[2], "matrix. Should be at least 5 x 2. Analyze more sets of SNPs and select more genomic features"),\"""" + output_path + """\",sep="\t", row.names=F, col.names=F) # Write clustering results            
            write.table(paste("ERROR: Cannot run correlation analysis on", dim(t5)[1], "x", dim(t5)[2], "matrix. Should be at least 5 x 2. Analyze more sets of SNPs and select more genomic features"),\"""" + pval_output_path + """\",sep="\t", row.names=F, col.names=F) # Write clustering results            

        }"""
    robjects.r(r_script)
    sys.stdout, sys.stderr = saved_stdout, saved_stderr
    # write matrices in json format     
    mat = open(output_path).read()      
    mat_pval = open(pval_output_path).read()        
    json_mat = []       
    json_mat.append({"log": False,"neg": "Underrepresented", "pos": "Overrepresented","name": "Pearsons","alpha":"","matrix": mat})      
    json_mat.append({"log": False,"neg": "Underrepresented", "pos": "Overrepresented","name": "P-value","alpha":"","matrix": mat_pval})      
    with open(".".join(output_path.split(".")[:-1]+ ["json"]),'wb') as f:       
        simplejson.dump(json_mat,f)
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

    # zip enrichment folder
    en_dir = os.path.join(outdir,'enrichment')
    z_en = zipfile.ZipFile(os.path.join(outdir,'enrichment.zip'),'a')
    for f in os.listdir(en_dir):
        z_en.write(os.path.join(en_dir,f),f)
    z_en.close()

    # zip annotation result folder if it exists
    anno_dir = os.path.join(outdir,'annotations')
    if os.path.exists(anno_dir):
        z_ano = zipfile.ZipFile(os.path.join(outdir,'annotations.zip'),'a')
        for f in os.listdir(anno_dir):
            z_ano.write(os.path.join(anno_dir,f),f)
        z_ano.close()

        

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
    output_files =  [os.path.join(outdir,x) for x in os.listdir(outdir) if x.endswith(".txt") or x.endswith(".pdf") or x.endswith('.zip')]
    fls = output_files
    for f in fls:
        tar.add(f,os.path.basename(f))
    tar.close()
    tar_file = open(tar_path,'rb')
    with gzip.open(tar_path+".gz","wb") as gz:
        gz.writelines(tar_file)
    tar_file.close()
    if os.path.exists(tar_path): os.remove(tar_path)

def validate_filenames(file_paths):
    ''' Checks if there are spaces before or after the file extension.
    EX. 'dir1/dir2/test .bed' is not valid. 'dir1/dir2/test.bed is valid.
    'dir1/dir2/test.bed .gz' is not valid.
    '''
    invalid = []
    for file in file_paths:
        for t in os.path.basename(file).split("."):
            if len(t.strip()) != len(t):
                # there are spaces before or after the '.'. Add the file to the list of invalids.
                invalid.append(os.path.basename(file))
    return invalid

def filter_score(gfs,pct_score):
    ''' Check the list of GFs to see if they exist in the pct_score filtered database.
    If it does not, then return the GF path to the non-filtered database. 
    If pct_score == '', return the GF paths from the non-filtered database.
    '''
    if pct_score == "":
        return gfs
    new_gfs = []
    for gf_path in gfs:
        gf_filtered_path = gf_path.replace('/grsnp_db/','/grsnp_db_{}/'.format(pct_score))
        if os.path.exists(gf_filtered_path):
            new_gfs.append(gf_filtered_path)
        else:
            new_gfs.append(gf_path)
    return new_gfs

def run_hypergeom(fois, gfs, bg_path,outdir,job_name="",zip_run_files=False,bkg_overlaps_path="",gr_data_dir = "" ,run_annotation=True,run_randomization_test=False,padjust="None",pct_score=""):
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

    logger.info("P_value adjustment used: {}".format(padjust))

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
        print fois
        logger.error(fois)
        gfs = [line for line in read_lines(gfs) if not line.endswith(".tbi")]
        # check if the GF exists in the database filtered by score
        gfs = filter_score(gfs,pct_score)
        # check if there are spaces in invalid parts of the file name
        invalid_names = validate_filenames(fois + gfs + [bg_path])
        if len(invalid_names) != 0:
            logger.error("The following file(s) have invalid file names, cannot have space before/in file extension:\n" + "\n".join(invalid_names))
            _write_progress("ERROR: Files have invalid filenames. See log file. Terminating run.")
            return         
        if bg_path.endswith(".tbi"):
            logger.error("Background has invalid extension (.tbi). Terminating run.")
            _write_progress("ERROR: Background has invalid extension (.tbi). Terminating run.")
            return

        foi_bg,good_fois = check_background_foi_overlap(bg_path,fois)   # Validate FOIs against background. Also get the size of the background (n_bgs)
        write_output("\t".join(map(base_name,good_fois))+"\n", matrix_outpath)
        write_output("\t".join(['foi_name', 'foi_obs', 'n_fois', 'bg_obs', 'n_bgs', 'odds_ratio', 'p_val','test_type','p_rand' if run_randomization_test else "",'p_mod' if run_randomization_test else ""]) + "\n",detailed_outpath)
        curprog,progmax = 0,len(gfs)
        # remove old detailed enrichment result files if they exit
        enr_path =  os.path.join(output_dir,"enrichment")
        for f in good_fois:
            f_path = os.path.join(enr_path, base_name(f)+'.txt')
            if os.path.exists(f_path): os.remove(f_path)
        _write_progress("Performing calculations on the background.")
        for gf in gfs: 
            current_gf = base_name(gf)      
            _write_progress("Performing Hypergeometric analysis for {}".format(base_name(gf)))
            str_pct_score = "Score threshold: NA"
            if '/grsnp_db_{}/'.format(pct_score) in gf:
                str_pct_score = "Score threshold: {}%".format(pct_score) 
            write_output("###"+base_name(gf)+"\t"+str_pct_score+"\t"+get_description(base_name(gf),trackdb)+"###"+"\n",detailed_outpath)
            res = get_overlap_statistics(gf,good_fois) 

            # calculate bg_obs
            bg_obs = get_bgobs(bg_path,gf,bkg_overlaps_path)
            if bg_obs == None: 
                logger.error("Skipping {}".format(gf))
                continue

            n_bgs = foi_bg[0]["indexregions"]  

            # run the enrichment analysis and output the matrix line for the current gf
            write_output("\t".join([base_name(gf)] + [str(p_value(res[i]["intersectregions"],res[i]["queryregions"],bg_obs,n_bgs ,good_fois[i],gf,bg_path,run_randomization_test)) for i in range(len(good_fois))])+"\n",matrix_outpath)
            curprog += 1
        # Adjust p values for the enrichment files
        list_enr =  [x  for x in os.listdir(enr_path)]

        for p in list_enr:
            adjust_detailed_pvalue(os.path.join(enr_path, p),padjust) 

        if len(gfs) > 1 and len(good_fois) > 1:
            # Adjust the pvalues
            adjust_pvalue(matrix_outpath,padjust)
            # Cluster the matrix
            clust_path =  cluster_matrix(matrix_outpath,os.path.join(outdir,"clustered.txt"))
            if len(gfs) > 4:               
                pearsons_cor_matrix(clust_path,outdir)
            else:
                json_mat = []
                json_mat.append({"log": False,"neg": "Underrepresented", "pos": "Overrepresented","name": "Pearsons","alpha":"","matrix": "ERROR:PCC matrix requires at least a 5 X 2 matrix."})      
                json_mat.append({"log": False,"neg": "Underrepresented", "pos": "Overrepresented","name": "P-value","alpha":"","matrix": "ERROR:PCC matrix requires at least a 5 X 2 matrix."})      
                with open(os.path.join(outdir,"pcc_matrix.json"),'wb') as f:   
                    simplejson.dump(json_mat,f)
                with open(os.path.join(os.path.join(outdir,"pcc_matrix.txt")),"wb") as wb:
                    wb.write("ERROR:PCC matrix requires at least a 5 X 2 matrix.")
        else:
            with open(os.path.join(outdir,"clustered.txt"),"wb") as wb:
                wb.write("ERROR:Clustered matrix requires at least a 2 X 2 matrix.")
        writer = open(os.path.join(outdir,"out.txt"),'wb')
        if run_annotation:
            annot_outdir = os.path.join(outdir,"annotations")
            if not os.path.exists(annot_outdir): os.mkdir(annot_outdir)
            curprog,progmax = 0,len(fois)
            for f in fois:                
                _write_progress("Running Annotation Analysis for {}.".format(base_name(f)))
                with open(os.path.join(annot_outdir,base_name(f) + ".txt"),"wb") as wr:
                    anot = get_annotation(f,gfs).split("\n")
                    anot[0] = anot[0].replace("Region\t\t","Region\t")
                    wr.write("\t".join(base_name(x) for x in anot[0].split("\t")) + "\tTotal")
                    for ind, a in enumerate(anot[1:]):
                        if a.strip() != "":
                            cur_row = a.split("\t")
                            wr.write("\n" + str(ind) + "|"+"\t".join(cur_row + [str(sum([int(x) for x in cur_row[1:] if x != ""]))]))                            
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
    valid_pv_adjust = ['bonferroni', 'holm', 'hochberg', 'hommel', 'BH', 'BY', 'fdr','None']
    parser = argparse.ArgumentParser(description="Enrichment analysis of several sets of SNPs (FOIs) files against several genomic features (GFs). Example: python hypergeom4.py foi_full_names.txt gf_full_names.txt /path_to_background/snp137.bed.gz")
    parser.add_argument("fois", nargs=1, help="Text file with paths to FOI files (unless -p used). Required") 
    parser.add_argument("gfs" ,nargs=1, help="Text file with pathrs to GF files (unless -p used). GF files may be gzipped. Required")
    parser.add_argument("bg_path", nargs=1, help="Path to background, or population of all SNPs. Required")
    parser.add_argument("--run_annotation" , "-a", help="Run annotation analysis", action="store_true" )
    parser.add_argument("--output_dir","-d", help="Directory to output the result to. Example: test_results. Default: current directory", default="")
    parser.add_argument("--pass_paths", "-p", help="Pass fois and gfs as comma separated paths. Paths are saved in .fois and .gfs file.", action="store_true")
    parser.add_argument("--pv_adjust", "-v",type=str, help="Which p-value adjustment method to use. Default: 'fdr'. Available (case-sensitive): "+', '.join(valid_pv_adjust), default="fdr")
        
    args = vars(parser.parse_args())  
    if args['pv_adjust'] not in valid_pv_adjust:
        print "ERROR: {} is not a valid p-value adjustment method.".format(args['pv_adjust'])
        sys.exit()

    if args["pass_paths"]: 
        gf = args["gfs"][0].split(",")      
        foi = args["fois"][0].split(",")        
        args["gfs"][0],args["fois"][0] = os.path.join(args["output_dir"],".gfs"),os.path.join(args["output_dir"],".fois")       
        with open(".gfs",'wb') as writer:       
            writer.write("\n".join(gf))     
        with open(".fois","wb") as writer:      
            writer.write("\n".join(foi))
    run_hypergeom(args["fois"][0],args["gfs"][0],args["bg_path"][0],args["output_dir"],"",False,"","",args["run_annotation"],run_randomization_test=False,padjust=args['pv_adjust'])

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


def _load_minmax(path):
    data = {}
    if not os.path.exists(path):
        return data
    score = [x for x in open(path).read().split("\n") if x != ""]
    for s in score:
        name,min_max = s.split('\t')
        data[name] = min_max
    return data