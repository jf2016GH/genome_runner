## Automates creating subgroups of genomic features into separate files/subfolders.
## Usage: Copy This Makefile and 'extract_UCSC.py' scripts into
## [databasedir]/custom_data/fois/hg19/, and run
##
## make
##
## It will create subfolders containing files with subgroups of genomic features,
## validate, sort, bgzip and tabix the files.
##
## Note: To create gwasCatalog subgroups, use https://github.com/mdozmorov/gwas2bed

all:	fois Cistrome.wget Epigenome.wget Motif.wget

fois:
		python extract_UCSC.py wgEncodeBroadHmmGm12878HMM name histoneMarks && \
		python extract_UCSC.py wgEncodeRegTfbsClusteredV3 name tfbsEncode && \
		python extract_UCSC.py gap type gapLocations && \
		python extract_UCSC.py coriellDelDup cellType coriellCNVs && \
		python extract_UCSC.py knownAlt name altSplicing && \
		python extract_UCSC.py wgRna type ncRnas && \
		python extract_UCSC.py tfbsConsSites name tfbsConserved && \
		python extract_UCSC.py dgvMerged varType genomicVariants && \
		python extract_UCSC.py nestedRepeats repClass repeats && \
		python extract_UCSC.py snp138Flagged func snpClinical && \
		echo "Data extraction complete"; \
		for file in `find . -type f -name "*.bed"`; do \
		f=`basename $$file`; d=`dirname $$file`; echo "Processing:" $$file; \
		cat $$file | grep "\bchr[0-9XYM][^_]\b" | \
		awk 'BEGIN {OFS="\t"} { if ( $$3 <= $$2) { print $$1, $$2, $$2+1, $$4, $$5, $$6 } \
		else { print $$1, $$2, $$3, $$4, $$5, $$6 } }' | \
		sort -k1,1 -k2,2n -k3,3n -u -o $$file && bgzip $$file && tabix $$file".gz"; \
		done

# Save sources of Cistrome web pages from http://cistrome.org/NR_Cistrome/ into corresponding HTML files
Cistrome.wget:	Cistrome.html
	grep ".bed" $< | tr '\"' '\t' | cut -f4 | sed 's/\.\///' | awk '{print "http://cistrome.org/NR_Cistrome/"$$0}' > $@ 

Epigenome.wget:	Epigenome.html
	grep ".bed" $< | tr '\"' '\t' | cut -f4 | sed 's/\.\///' | awk '{print "http://cistrome.org/NR_Cistrome/"$$0}' > $@ 

Motif.wget:	Motif.html
	grep ".bed" $< | tr '\"' '\t' | cut -f4 | sed 's/\.\///' | awk '{print "http://cistrome.org/NR_Cistrome/"$$0}' > $@ 

