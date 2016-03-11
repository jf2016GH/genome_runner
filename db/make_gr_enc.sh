#!/bin/bash

# ENCODE enrichment analysis

# Define variables
# Path to BED files to test
DIR=/home/mdozmorov/work/gwas2bed/gwasCatalog/
# Path to GenomeRunner's database
DB=/home/mdozmorov/db_5.00_07-22-2015
# Path to results folder
RES=/home/mdozmorov/work/gwas2bed/gwasCatalog/data.gr
# Path to the background file
BKG=/home/mdozmorov/work/gwas2bed/gwasCatalog/all_current_gwas_catalog.bed.gz

# Make the results folder
mkdir -p $RES

# Create text files with FOIs
find $DIR/ -type f -name "*.txt" | sort > fois.txt

# Enrichment analysis. GF files should be created
for file in `ls /home/mdozmorov/work/gfs* | grep enc | sort`; do 
	clear;
	echo $file;
	GF=$file;  
	DIROUT=$RES"/gr_"`basename $GF .txt | sed 's/gfs_//'`;
	gr -g hg19 -d $DB -r $DIROUT fois.txt $GF $BKG;
	sleep 1;
done

