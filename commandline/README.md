`make_bkg_hg19.sh` - Download snp142 and snp142Common files from the UCSC database, and process them into BED format. Non-canonical chromosomes are removed.

`make_download_LD.sh` - Download LD scores for EUR population from [https://data.broadinstitute.org/srlab/BEAGLE/1kG-beagle-release3/pairwise_ld/r2_ge_0.8/](https://data.broadinstitute.org/srlab/BEAGLE/1kG-beagle-release3/pairwise_ld/r2_ge_0.8/)

# How to calculate % of genome coverage for a BED file

1. Get some BED files:

		wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_GM12878_subcompartments.bed.gz
		zcat < GSE63525_GM12878_subcompartments.bed.gz | cut -f1-4 | grep -v NA | awk 'BEGIN {OFS="\t"} {print $0 >> ("TAD_"$4".bed"); close("TAD_"$4".bed")}'
		for file in TAD_*; do bedtools sort -i $file > tmp.bed && mv tmp.bed $file; done

2. Get genomic coordinates of autosomal chromosomes

		mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg19.chromInfo" | awk '{OFS="\t"} {print $1,"0",$2}' | sed '1d' | grep -v "_" | bedtools sort -i - > hg19.bed

3. Calculate coverage by chromosome (last column contains fraction of the first file overlapping the second (chromosome) file)

		coverageBed -a TAD_B4.bed -b hg19.bed

or, the total percent

		coverageBed -a TAD_B4.bed -b hg19.bed | cut -f7 | datamash sum 1