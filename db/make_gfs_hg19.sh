#!/bin/bash

# Create a list of GFs
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/hg19/genes/cytobands/ -type f -name "*.bed.gz" | sort > gfs_cytobands.txt
find /home/mdozmorov/db_5.00_07-22-2015/custom_data/gfs/hg19/GERP/ -type f -name "*.bed.gz" | sort > gfs_GERP.txt
find /home/mdozmorov/db_5.00_07-22-2015/custom_data/gfs/hg19/CpG/ -type f -name "*.bed.gz" | sort > gfs_CpG.txt
find /home/mdozmorov/db_5.00_07-22-2015/custom_data/gfs/hg19/gwasCatalog/ -type f -name "*.bed.gz" | sort > gfs_gwasCatalog.txt
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/hg19/genes/long_noncoding/ -type f -name "*.bed.gz" | sort >> gfs_genes.txt
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/hg19/genes/protein_coding/ -type f -name "*.bed.gz" | sort >> gfs_genes.txt
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/hg19/genes/pseudogene/ -type f -name "*.bed.gz" | sort >> gfs_genes.txt
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/hg19/genes/short_noncoding/ -type f -name "*.bed.gz" | sort >> gfs_genes.txt
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/hg19/genes/KEGG/ -type f -name "*.bed.gz" | sort >> gfs_KEGG.txt
find /home/mdozmorov/db_5.00_07-22-2015/custom_data/gfs/hg19/DGV/ -type f -name "*.bed.gz" | sort > gfs_DGV.txt
find /home/mdozmorov/db_5.00_07-22-2015/custom_data/gfs/hg19/H3K4me3/ -type f -name "*.bed.gz" | sort > gfs_H3K4me3.txt
find /home/mdozmorov/db_5.00_07-22-2015/custom_data/gfs/hg19/nestedRepeats/ -type f -name "*.bed.gz" | sort > gfs_nestedRepeats.txt
find /home/mdozmorov/db_5.00_07-22-2015/custom_data/gfs/hg19/VMRs/ -type f -name "*.bed.gz" | sort > gfs_VMRs.txt
find /home/mdozmorov/db_5.00_07-22-2015/custom_data/gfs/hg19/super_enhancers/ -type f -name "*.bed" | sort > gfs_super_enhancers.txt
find /home/mdozmorov/db_5.00_07-22-2015/custom_data/gfs/hg19/UCNEs/ -type f -name "*.bed.gz" | sort > gfs_UCNEs.txt
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/hg19/FANTOM/ -type f -name "*.bed.gz" | sort >> gfs_FANTOM.txt

# ENCODE GFs
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/hg19/ENCODE/TFBS_clustered/ -type f -name "*.bed.gz" | sort > gfs_encTFBS_clustered.txt
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/hg19/ENCODE/TFBS_cellspecific/ -type f -name "*.bed.gz" | sort > gfs_encTFBS_cellspecific.txt
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/hg19/ENCODE/Histone/ -type f -name "*.bed.gz" | sort > gfs_encHistone.txt
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/hg19/ENCODE/DNase/ -type f -name "*.bed.gz" | sort > gfs_encDNase.txt
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/hg19/ENCODE/chromStates/ -type f -name "*.bed.gz" | sort > gfs_encchromStates.txt

# Roadmap GFs
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/hg19/ROADMAP/chromStates15/ -type f -name "*.bed.gz" | sort > gfs_rdmchromStates15.txt
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/hg19/ROADMAP/chromStates18/ -type f -name "*.bed.gz" | sort > gfs_rdmchromStates18.txt
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/hg19/ROADMAP/chromStates25/ -type f -name "*.bed.gz" | sort > gfs_rdmchromStates25.txt
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/hg19/ROADMAP/DNase_bPk-processed/ -type f -name "*.bed.gz" | sort > gfs_rdmDNase_bPk-processed.txt
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/hg19/ROADMAP/DNase_gPk-imputed/ -type f -name "*.bed.gz" | sort > gfs_rdmDNase_gPk-imputed.txt
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/hg19/ROADMAP/DNase_nPk-imputed/ -type f -name "*.bed.gz" | sort > gfs_rdmDNase_nPk-imputed.txt
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/hg19/ROADMAP/DNase_nPk-processed/ -type f -name "*.bed.gz" | sort > gfs_rdmDNase_nPk-processed.txt
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/hg19/ROADMAP/Histone_bPk-processed/ -type f -name "*.bed.gz" | sort > gfs_rdmHistone_bPk-processed.txt
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/hg19/ROADMAP/Histone_gPk-imputed/ -type f -name "*.bed.gz" | sort > gfs_rdmHistone_gPk-imputed.txt
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/hg19/ROADMAP/Histone_gPk-processed/ -type f -name "*.bed.gz" | sort > gfs_rdmHistone_gPk-processed.txt
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/hg19/ROADMAP/Histone_nPk-imputed/ -type f -name "*.bed.gz" | sort > gfs_rdmHistone_nPk-imputed.txt
find /home/mdozmorov/db_5.00_07-22-2015/grsnp_db/hg19/ROADMAP/Histone_nPk-processed/ -type f -name "*.bed.gz" | sort > gfs_rdmHistone_nPk-processed.txt

