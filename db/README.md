Misc. scripts and files for working with database of genome annotation files
================================================

- `CISTROME`
	- `Makefile.cistrome`, `Cistrome.*`, `Epigenome.*`, `Motif.*` - Makefile and prerequisites to download and create Cistrome database for GenomeRunner

- `autorsync.sh`, `autorsync1.sh` - downloads whole UCSC database using rsync, restarting rsync if it breaks

- `blacklisted.txt` - genome annotation features that should be ignored

- `extract_UCSC.py` - split multiple-feature genome annotation supertracks into separate feature-specific features

- `Makefile` - Makefile for creating hg19 and mm9 custom genomic features. Currently supported: CpG BED files (`make cpg_hg19` and `make cpg_mm9`), from [http://rafalab.jhsph.edu/CGI/](http://rafalab.jhsph.edu/CGI/); Cytobands (`make cytoband_hg19`), Genomic Evolutionary Rate Profiling (GERP) elements (`make GERP_hg19` and `make GERP_mm9`)

- `Makefile.parallel`, `make-matrix.py` - using an example of parallel processing from [BedTools tutorial](https://github.com/arq5x/tutorials/blob/master/bedtools.md) to get Jaccard statistics for overlaps among multiple BED files.