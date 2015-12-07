all:	cpg_hg19 cpg_mm9 cytoband_hg19

# Create CpG BED files, from http://rafalab.jhsph.edu/CGI/

model-based-cpg-islands-hg19.txt:
	wget http://rafalab.jhsph.edu/CGI/model-based-cpg-islands-hg19.txt

cpg_hg19:	model-based-cpg-islands-hg19.txt
	sed '1d' $< | grep "\bchr[0-9XYM][^_]\b" | awk 'BEGIN {OFS="\t"} { if ( $$3 <= $$2) { print $$1, $$2, $$2+1 } else { print $$1, $$2, $$3 } }' | sort -k1,1 -k2,2n -k3,3n | uniq > CpG-$@.bed;
	bgzip CpG-$@.bed;
	tabix CpG-$@.bed.gz


model-based-cpg-islands-mm9.txt:
	wget http://rafalab.jhsph.edu/CGI/model-based-cpg-islands-mm9.txt

cpg_mm9:	model-based-cpg-islands-mm9.txt
	sed '1d' $< | grep "\bchr[0-9XYM][^_]\b" | awk 'BEGIN {OFS="\t"} { if ( $$3 <= $$2) { print $$1, $$2, $$2+1 } else { print $$1, $$2, $$3 } }' | sort -k1,1 -k2,2n -k3,3n | uniq > CpG-$@.bed;
	bgzip CpG-$@.bed;
	tabix CpG-$@.bed.gz

cytoBand.txt.gz:
	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz

cytoband_hg19:	cytoBand.txt.gz
	gzip -d $<;
	while read p; do echo $$p | sed "s/\./\_/g" | awk 'BEGIN {OFS="\t"} {print $$1,$$2,$$3,$$4,$$5 >> ($$1"_"$$4".bed"); close($$1"_"$$4".bed")}' ; done <cytoBand.txt;
	rm cytoBand.txt;
	for file in *.bed; do bgzip $$file; tabix $$file.gz; done
