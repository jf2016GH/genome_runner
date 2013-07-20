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


Interval = collections.namedtuple("Interval", "chrom,start,end")

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
current_gf,current_foi = "", "" # used for annotation analysis
#  There is one entry for each FOI.bed. "FOIname": FOI_Annotation
annotations = {}

def read_intervals(path, background=None, snp_only=False, d_path=None):
    with (gzip.open(path,"rb") if path.endswith(".gz") else open(path,"rb")) as h:
        intervals = []
        num_before, num_zero,num_great,num_back,num_chrom = 0, 0, 0, 0,0

        for i,line in enumerate(h):
            num_before += 1
            fields = line.strip().split("\t")
#            print fields
            chrom, start, end = fields[:3]
            start = int(start)
            end = int(end)

            # exclude SNPs from weired chromosomes
            if "_" not in chrom:
                # convert regions to SNPs
                if snp_only and ((end - start) > 1):
                    logger.warning("\t"+" ".join([chrom,str(start),str(end)]) + " in\t{}\tis not a SNP. Shortening feature to be SNP".format(path))
                    start = (end-start)/2
                    end, num_great = start + 1, num_great + 1

                # convert '0 length' regions to SNPs
                if snp_only and ((end - start) == 0):
                    logger.warning("\t"+" ".join([chrom,str(start),str(end)]) + " in\t{}\tis zero length. Converting to SNP".format(path))
                    end, num_zero = start + 1, num_zero + 1

                interval = Interval(chrom, int(start), int(end))
                # output warning for features not in background
                if background and not background.query(interval):
                    logger.warning("\t"+" ".join([chrom,str(start),str(end)]) + " in\t{}\tdoes not overlap with the background. Continuing...".format(path))                
                    num_back += 1
                intervals.append(interval)
            else:
                logger.warning("\t"+" ".join([chrom,str(start),str(end)]) + " in\t{}\tis on strange chromosome, excluding".format(path))
                num_chrom += 1

        # output grooming summary
        if d_path:
            str_groom_sum = """\n{}
            {} SNPs before grooming.
            {} SNPs with length of 0 converted to SNP
            {} SNPs with length > 1 converted to SNP
            {} SNPs on strange chromosomes and excluded
            {} SNPs after grooming.""".format(os.path.split(path)[-1], num_before,num_zero,num_great,num_chrom,len(intervals))
            _write_head(textwrap.dedent(str_groom_sum),d_path)
    return intervals

class IntervalSet(object):
    def __init__(self, intervals, snp_only=False):
        self._trees = {}
        self._intervals = list(intervals)
        for i,iv in enumerate(intervals):
            if not iv.chrom in self._trees:
                self._trees[iv.chrom] = IntervalTree()
            self._trees[iv.chrom].add(iv.start, iv.end, i)

    def query(self, interval):
        '''Checks for overlap of single SNP.
        '''
        t = self._trees.get(interval.chrom)
        # Found a chromosome
        if t: 
            return t.find(interval.start, interval.end)  

    def query_annotation(self, interval):
        '''Checks for overlap of single SNP. Does Annotation
        '''
        t = self._trees.get(interval.chrom)
        if t: 
            found = t.find(interval.start, interval.end)   
            # relies on a global dictionary to annotation results
            global annotations
            annotations[current_foi].add_hit(current_gf,int(found != []))
            return found
        # record a miss for annotation in case no GF is on the same chrom as the FOI
        else:
            global annotations
            annotations[current_foi].add_hit(current_gf,int(0))


    def n_overlaps(self, intervals,annotate=False):
        '''Gets the total number over laps. Records annotation results if annotation=True
        '''
        if annotate:
            return sum([1 if self.query_annotation(iv) else 0 for i,iv in enumerate(intervals)])
        else:
            return sum([1 if self.query(iv) else 0 for iv in intervals])

    def __len__(self):
        return len(self._intervals)

def p_value(gf, fois, bgs, foi_name,annotate):
    "Return the signed log10 p-value of intersection with the given interval set."
    global current_foi
    current_foi = foi_name
    if not foi_name in annotations.keys(): annotations[foi_name] = FOI_Annotation(foi_name,fois)  # create new annotation entry
    foi_obs = gf.n_overlaps(fois,annotate) # number of FOIs overlapping with a GF
    bg_obs = gf.n_overlaps(bgs) # number of spot bkg overlapping with a GF
    n_fois, n_bgs = len(fois), len(bgs)
    ctable = [[foi_obs, n_fois-foi_obs],
              [bg_obs-foi_obs,n_bgs-n_fois-(bg_obs-foi_obs)]]
    odds_ratio, pval = scipy.stats.fisher_exact(ctable)
    sign = 1 if (odds_ratio < 1) else -1
    write_output("\t".join(map(str, [foi_name.rpartition('/')[-1], foi_obs, n_fois, bg_obs, n_bgs, "%.2f" % odds_ratio, "%.2f" % pval])) + "\n",detailed_outpath)
    try:        
        return sign * math.log10(pval)        
    except ValueError as e:
        logger.error(e)
        return sign * sys.float_info.min_10_exp



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
                    h = heatmap.2(pt5,  hclustfun=function(m) hclust(m,method="average"),  distfun=function(x) dist(x,method="euclidean"), cexCol=1, cexRow=1)
                    write.table(t(h$carpet),"{}",sep="\t")""".format(input_path,output_path)
    robjects.r(r_script)
    return output_path    

def pearsons_cor_matrix(matrix_path,out_dir):
    output_path = os.path.join(out_dir,".cor")
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
                    write.table(t(h$carpet),"{}",sep="\t")""".format(matrix_path,output_path)
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
                write.table(pt5.prob,\""""+output_path+".pvalue"+"""\",sep="\t")"""
    #robjects.r(r_script)
    r_script = """t5 = read.table(\""""+matrix_path+"""\") 
        library("Hmisc")
        p5<-rcorr(as.matrix(t5))
        p5[[1]]
        p5[[3]]
        h<-heatmap.2(as.matrix(p5[[1]]))
        write.table(h$carpet,\"""" + output_path + """\",sep="\t")
        write.table(p5[[3]][h$rowInd, h$colInd],\""""+output_path+".pvalue"+ """\",sep="\t")""" 
    robjects.r(r_script)
    return output_path   



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
    return os.path.basename(path).split(".")[0]

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


def run_hypergeom(fois, gfs, bg_path,outdir,job_name="",zip_run_files=False,run_annotation=False):
    sett_path = os.path.join(outdir,".settings")
    logger_path = os.path.join(outdir,'.log')
    trackdb = []
    if os.path.exists(sett_path):        
        with open(sett_path) as re:
            organism = [x.split("\t")[1] for x in re.read().split("\n") if x.split("\t")[0] == "Organism:"][0]
            trackdb = bedfilecreator.load_tabledata_dumpfiles(os.path.join("data",organism,"trackDb"))
    # set output settings
    global detailed_outpath,matrix_outpath, progress_outpath, curprog, progmax,current_gf
    detailed_outpath =  os.path.join(outdir, "detailed.gr")
    matrix_outpath = os.path.join(outdir,"matrix.gr")
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
    bg = read_intervals(bg_path)
    bg_iset = IntervalSet(bg)

    _write_head("\n\n#Detailed log report#\n",logger_path)
    foi_sets = dict((path,read_intervals(path, snp_only=True, background=bg_iset,d_path=logger_path)) for path in fois)
    _write_head("#Grooming Summary#",logger_path)

    write_output("\t".join(map(base_name,fois))+"\n", matrix_outpath)
    write_output("\t".join(['foi_name', 'foi_obs', 'n_fois', 'bg_obs', 'n_bgs', 'odds_ratio', 'p_val']) + "\n",detailed_outpath)
    curprog,progmax = 0,len(gfs) 
    try:
        for gf in gfs:        
            current_gf = base_name(gf)
            curprog += 1
            _write_progress("Performing Hypergeometric analysis for {}".format(base_name(gf)))
            gf_iset = IntervalSet(read_intervals(gf))            
            write_output("###"+base_name(gf)+"\t"+get_description(base_name(gf),trackdb)+"###"+"\n",detailed_outpath)
            write_output("\t".join([base_name(gf)] + [str(p_value(gf_iset, foi_sets[foi], bg, foi,run_annotation)) for foi in fois])+"\n",matrix_outpath)
        if len(gfs) > 1 and len(fois) > 1:
            cluster_matrix(matrix_outpath,os.path.join(outdir,"matrix_clustered.gr"))
            pearsons_cor_matrix(os.path.join(outdir,"matrix_clustered.gr"),outdir)
        else:
            with open(os.path.join(outdir,"matrix_clustered.gr"),"wb") as wb:
                wb.write("Clustered matrix requires at least a 2 X 2 matrix.")
        if run_annotation:
            _write_progress("Outputting annotation data")
            with open(os.path.join(outdir, "annotations.gr"),"w") as a_out:
                for k,v in annotations.iteritems():
                    a_out.write(v.return_str_matrix() +"\n")
        _write_progress("Preparing run files for download")
        _zip_run_files(fois,gfs,bg_path,outdir,job_name)
        _write_progress("Analysis Completed")       
    except Exception, e: 
        logger.error( traceback.print_exc())
        _write_progress("Run crashed. See end of log for details.")

def get_description(gf,trackdb):
    desc = [x["longLabel"] for x in trackdb if x["tableName"] == gf]
    if len(desc) is not 0: return desc[0]
    else: return "No Description"




def _zip_run_files(fois,gfs,bg_path,outdir,job_name=""):
    '''
    File paths of FOIs and GFs as a list. Gathers all the files together in one zipped file
    '''    
    f = open(os.path.join(outdir,".log"))
    f_log = f.read()
    f.close()
    path_settings,f_sett =os.path.join(outdir,".settings"),""
    if os.path.exists(path_settings):
        f = open(path_settings)
        f_sett = f.read() + "\n###LOG###\n"
        f.close()
    new_log_path = os.path.join(outdir,".details")
    new_log = open(new_log_path,'wb')
    new_log.write(f_sett+f_log)
    new_log.close()
    tar_path = os.path.join(outdir,'GR_Runfiles_{}.tar'.format(job_name))
    tar = tarfile.TarFile(tar_path,"a")    
    output_files =  [os.path.join(outdir,x) for x in os.listdir(outdir) if x.endswith(".gr")]
    fls = output_files + [new_log_path]
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
    parser.add_argument("fois", nargs=1, help="Text file with FOI file names (SNPs only).") 
    parser.add_argument("gfs", nargs=1, help="Text file with GF file, gzipped.") 
    parser.add_argument("bg_path", nargs=1, help="Path to spot background file (SNPs only).")
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
            wb.write("Clustered matrix requires at least a 2 X 2 matrix.")


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


# 'gf': name of the GF. 
# Each entry in 'hits' represents a single SNP.
# The order must match the foi.bed.
Annotation = collections.namedtuple("Annotation", "gf,hits") 

class FOI_Annotation:
    '''Stores the name of the FOI.bed and creates a 
    list of 'Annotation' for each GF run.
    'Annotation' stores GF.bed name and a list of booleans, each boolean
     represents whether the corresponding SNP overlapped with the GF.
    '''
    def __init__(self,foiname,foi_intervals):
        self._foiname = foiname
        self._fois = foi_intervals
        self._results = [] # stores results as Annotation entries
    
    def add_hit(self, gf,hit):
        ''' 'gf' is the name of the GF.
            'hit' is boolean of whether SNP overlapped.  
            Appends 'hit' to existing 'Annotation' if one exists
            for the 'gf'. Otherwise, creates new 'Annotation'          
        '''
        ind = [i for i,x in enumerate(self._results) if x.gf==gf]

        if len(ind) != 0:
            self._results[ind[0]].hits.append(hit)
        else:
            self._results.append(Annotation(gf,[hit]))

    def return_str_matrix(self):
        '''Returns a string of the formated annotation results matrix.
        '''
        res = "###{}###\n".format(self._foiname) 
        res += "FOI\t" + "\t".join([x.gf for x in self._results]) + "\tTotal\n"
        for i,v in enumerate(self._fois):         
            res += "|".join(str(x) for x in v) + "\t"
            # print out each hit for each GF record for the current SNP
            hits = [x.hits[i] for x in self._results ]
            res += "\t".join(str(x) for x in hits) + "\t" + str(sum(hits)) + "\n"
        return res
