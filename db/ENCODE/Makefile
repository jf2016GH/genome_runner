#TYPE = Hist
TYPE = Tfbs

all:	stats_ENCODE_$(TYPE).txt

SHELL = /bin/bash
DB = /home/genomerunner/db_3.00_12-01-2014
ORG = hg19

stats_ENCODE_$(TYPE).txt:
	mkdir 111; \
	mkdir 222; \
	for file1 in `find $(DB)/grsnp_db/$(ORG)/ENCODE/ -type f -name "*.bed.gz" | grep -i $(TYPE) |sort`; do \
		for file2 in `find $(DB)/grsnp_db/$(ORG)/ENCODE/ -type f -name "*.bed.gz" | grep -i $(TYPE) | sort`; do \
			f1bed=`basename $$file1`; \
			f1bed=$${f1bed%.gz}; \
			f2bed=`basename $$file2`; \
			f2bed=$${f2bed%.gz}; \
			zcat $$file1 > 111/$$f1bed ; \
			zcat $$file2 > 222/$$f2bed && \
			echo -n $$f1bed"	"$$f2bed"	" >> $@; \
			bedtools jaccard -a 111/$$f1bed -b 222/$$f2bed | sed '1d' >> $@ && \
			echo "Processing " $$f1bed "   " $$f2bed ; \
			rm 111/$$f1bed; \
			rm 222/$$f2bed; \
		done; \
	done; \
	rm -r 111; \
	rm -r 222


