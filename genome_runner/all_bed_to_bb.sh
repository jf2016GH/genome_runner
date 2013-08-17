for file in `find . -type f -name '*gz'`;
do 
d=`dirname $file`;
f=`basename $file`;
f=${f%.gz};
gunzip $file && ~/bedToBigBed $d"/"$f ~/hg19.chromInfo.txt $d"/"$f".bb" && rm $d"/"$f;
done 
