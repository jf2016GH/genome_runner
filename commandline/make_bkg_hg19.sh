
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/snp142.txt.gz # 2.4G
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/snp142Common.txt.gz # 522M

FILE=snp142.txt.gz

zcat $FILE | grep "\bchr[0-9XYM][^_]\b" | awk 'BEGIN {OFS="\t"} { if ( $3 <= $2) { print $1, $2, $2+1 } else { print $1, $2, $3 } }' | sort -k1,1 -k2,2n -k3,3n | uniq

zcat $FILE | grep "\bchr[0-9XYM][^_]\b" | awk 'BEGIN {OFS="\t"} { if ( $3 <= $2) { print $1, $2, $2+1, $4 } else { print $1, $2, $3, $4 } }' | sort -k4,4 | uniq $FILE".bed"

zcat snp142Common.txt.gz | awk 'BEGIN {OFS="\t"} { { if ( $4 <= $3 ) { print $2, $3, $3+1, $5 } else { print $2, $3, $4, $5 } }' | grep "\bchr[0-9XYM][^_]\b" | sort -k4,4 | uniq > snp142Common.bed
