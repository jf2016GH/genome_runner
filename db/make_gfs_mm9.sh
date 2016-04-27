#!/bin/bash

# Create a list of GFs
find /home/mdozmorov/db_5.00_07-22-2015/custom_data/gfs/mm9/GERP/ -type f -name "*.bed.gz" | sort > gfs_GERP_mm9.txt
find /home/mdozmorov/db_5.00_07-22-2015/custom_data/gfs/mm9/CpG/ -type f -name "*.bed.gz" | sort > gfs_CpG_mm9.txt

# ENCODE GFs
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/mm9/ENCODE/DNase/ -type f -name "*.bed.gz" | sort > gfs_encDNase_mm9.txt
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/mm9/ENCODE/Histone/ -type f -name "*.bed.gz" | sort > gfs_encHistone_mm9.txt
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/mm9/ENCODE/TFBS/ -type f -name "*.bed.gz" | sort > gfs_encTFBS_mm9.txt

# CISTROME GFs
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/mm9/CISTROME/Cistrome/ -type f -name "*.bed.gz" | sort > gfs_cisCistrome_mm9.txt
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/mm9/CISTROME/Epigenome/ -type f -name "*.bed.gz" | sort > gfs_cisEpigenome_mm9.txt
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/mm9/CISTROME/Motif/ -type f -name "*.bed.gz" | sort > gfs_cisMotif_mm9.txt

