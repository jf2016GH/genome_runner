Misc. scripts and files for working with database of genome annotation files
================================================

- `CISTROME`
	- `Makefile.cistrome`, `Cistrome.*`, `Epigenome.*`, `Motif.*` - Makefile and prerequisites to download and create Cistrome database for GenomeRunner

- `autorsync.sh`, `autorsync1.sh` - downloads whole UCSC database using rsync, restarting rsync if it breaks

- `blacklisted.txt` - genome annotation features that should be ignored

- `extract_UCSC.py` - split multiple-feature genome annotation supertracks into separate feature-specific features

- `Makefile` - Makefile for creating hg19 and some mm9 custom genomic features. Currently supported: 

- CpG BED files (`make cpg_hg19` and `make cpg_mm9`), from [http://rafalab.jhsph.edu/CGI/](http://rafalab.jhsph.edu/CGI/) (link is broken);
- Cytobands (`make cytoband_hg19`);
- Genomic Evolutionary Rate Profiling (GERP) elements (`make GERP_hg19` and `make GERP_mm9`);
- CpGs and Variably Methylated Regions (VMRs) across 54 normal human cell types (See Gu J, et. al. ["Mapping of Variable DNA Methylation across Multiple Cell Types Defines a Dynamic Regulatory Landscape of the Human Genome."](https://www.ncbi.nlm.nih.gov/pubmed/26888867) G3 (Bethesda) 2016 for more details, (`make VMC_hg19` and `make VMR_hg19`;
- UltraConserved Noncoding Elements (UCNEs from the [UCNE base](http://ccg.vital-it.ch/UCNEbase/)). See Dimitrieva, S. and Bucher, P. ["UCNEbase – a database of ultra-conserved non-coding elements and genomic regulatory blocks."](https://www.ncbi.nlm.nih.gov/pubmed/23193254), 2013 for more details. (`make UCNE.bed`).



- `Makefile.parallel`, `make-matrix.py` - using an example of parallel processing from [BedTools tutorial](https://github.com/arq5x/tutorials/blob/master/bedtools.md) to get Jaccard statistics for overlaps among multiple BED files.make_gfs.sh

- `make_gfs.sh` - script to make GF files, containing lists of categories of GFs

- `make_gr_misc.sh` - script to run GR analysis on miscellaneous categories of GFs and custom data

- `make_gr_enc.sh` - script to run GR analysis on the ENCODE categories of GFs

- `make_gr_rdm.sh` - script to run GR analysis on the Roadmap categories of GFs

