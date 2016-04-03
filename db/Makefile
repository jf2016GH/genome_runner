# all:	cpg_hg19 cpg_mm9 cytoband_hg19
# all:	GERP_hg19 GERP_mm9
# all:	VMC_hg19 VMR_hg19
all:	UCNE.bed

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

# Generate GERP conserved elements from http://mendel.stanford.edu/SidowLab/downloads/gerp/

hg19.GERP_elements.tar.gz:
	wget http://mendel.stanford.edu/SidowLab/downloads/gerp/hg19.GERP_elements.tar.gz

GERP_hg19:	hg19.GERP_elements.tar.gz
	tar -zxvf $<;
	for file in hg19*.txt; do chr=`echo $$file | sed 's/hg19_//' | sed 's/_elems.txt//'`; awk -v chr=$$chr '{FS="\t"} {OFS="\t"} {print chr,$$1,$$2,$$3,$$4}' $$file; done | bedtools sort -i  > $@.bed;
	bgzip $@.bed;
	tabix $@.bed.gz

mm9.GERP_elements.tar.gz:
	wget http://mendel.stanford.edu/SidowLab/downloads/gerp/mm9.GERP_elements.tar.gz

GERP_mm9:	mm9.GERP_elements.tar.gz
	tar -zxvf $<;
	for file in mm9*.txt; do chr=`echo $$file | sed 's/mm9_//' | sed 's/_elems.txt//'`; awk -v chr=$$chr '{FS="\t"} {OFS="\t"} {print chr,$$1,$$2,$$3,$$4}' $$file; done | bedtools sort -i  > $@.bed;
	bgzip $@.bed;
	tabix $@.bed.gz

# Variably methylated regions
vmc.txt.gz:
	wget http://epigenome.wustl.edu/methylomes/download/vmc.txt.gz

VMC_hg19:	vmc.txt.gz
	zcat $< | bedtools sort -i - > $@.bed;
	bgzip $@.bed;
	tabix $@.bed.gz

# Variably methylated CpG sites
vmr.txt.gz:
	wget http://epigenome.wustl.edu/methylomes/download/vmr.txt.gz

VMR_hg19:	vmr.txt.gz
	zcat $< | cut -f1-3 | bedtools sort -i - > $@.bed;
	bgzip $@.bed;
	tabix $@.bed.gz

# UltraConserved Noncoding elements
UCNE.bed:
	wget http://ccg.vital-it.ch/UCNEbase/data/download/ucnes/hg19_UCNE_coord.bed
	bedtools sort -i hg19_UCNE_coord.bed > $@
# wget http://ccg.vital-it.ch/UCNEbase/data/download/ucnes/UCNE_paralog_regions.txt
# sed '1d' UCNE_paralog_regions.txt | cut -f3-5 | grep -v "_" | bedtools sort -i - | uniq > UCNE_paralog_regions.bed


