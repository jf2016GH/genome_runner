import textwrap,os


feature_data = {
'foi_overlap_test': '''chr1	886763	886763	(blank): A single point FOI (Feature of Interest) that lies outside of the GF(Genomic Feature) region
chr1	886783	886784	(SNPWithin): a SNP that lies inside of the GF region
chr1	886798	886799	(SNPWithin): a SNP that is equal to the right endpoint of the GF
chr1	886773	886774	(SNPWithin): a SNP that is equal to the left endpoint of the GF
chr1	886773	886798	(Exact): the startpoint and endpoint are equal to the startpoint and endpoint of the GF
chr1	886763	886790	(Overlap): the startpoint of the FOI is less than the GF startpoint and the FOI endpoint falls between the startpoint and the endpoint of the GF region
chr1	886783	886808	(Overlap): the startpoint of the FOI is between the endpoint and startpoint of the GF and the endpoint of the FOI is greater than the endpoint of the GF
chr1	886763	886773	"(Overlap, value = 1bp)  the start of FOI is less than the start of the GF and the endpoint of the FOI is equal to the startpoint of the GF region"
chr1	886798	886808	"(Overlap, value = 1bp) the start of FOI is equal to the endpoint of GF and the end of FOI is greater than the endpoint of the GF region"
chr1	886763	886808	(Overhang): the startpoint of the FOI is less than the startpoint of the GF region and the endpoint of the FOI is greater than the endpoint of the GF region
chr1	886783	886788	(Within): the startpoint of the FOI is greater than the start of the GF and the endpoint of the FOI is less than the endpoint of the GF
chr1	886763	886769	(blank): occurs when the startpoint and endpoint of the FOI is less than the startpoint of the GF
chr1	886808	886818	(blank): occurs when the startpoint and endpoint of the FOI is greater than the endpoint of the GF
chr1	886773	886788	(Within) the startpoint of the FOI is equal to the start of the GF and the endpoint of the FOI is less than the endpoint of the GF
chr1	886783	886798	(Within) the startpoint of the FOI is greater than the start of the GF and the endpoint of the FOI is equal to the endpoint of the GF
chr1	886773	886808	(Overhang): the startpoint of the FOI is equal to the startpoint of the GF and the endpoint of the FOI is greater then the endpoint of the FOI
chr1	886763	886798	(Overhang): the startpoint of the FOI is less than the startpoint of the GF and the endpoint of the FOI is equal to the endpoint of the GF  
chr1	1109337	1182285	"(Overlap,Overhang,Overlap): here the startpoint of the FOI is equal to the endpoint of  a GF region to the left (Overlap, value = 0), the FOI region completely overlaps a middle GF region (overhang), and overlaps a region on the right (Overlap)"
chr1	1109325	1182277	"(Overlap, Overhang, Overlap): here the startpoint of the FOI is in the middle of a GF region on the left resulting in an Overlap, the FOI completely overlaps a middle GF region (overhang), and the endpoint of the FOI is equal to the startpoint of a third region on the right (counted as blank)"
chr1	1109337	1182277	"(Overlap,Overhang,Overlap): here the startpoint of the FOI is equal to the endpoint of a GF region on the right, the FOI completely overlaps a middle GF region(Overhang), and the endpoint of the FOI is equal to the startpoint of the GF region on the right (Overlap, value =0)"
chr1	1109325	1158465	"(Overlap,Overlap): The FOI region overlaps two different GF regions"
chr1	1158465	1182285	"(Overlap,Overlap): The FOI region overlaps two different GF regions"
''',

'gf_overlap_test': '''chr1	886773	886798	608_0_-_96	96	-
chr1	888417	888435	617_0_+_156	156	+
chr1	1094265	1094313	1087_0_+_62	62	+
chr1	1109319	1109337	1109_0_+_111	111	+
chr1	1158458	1158482	1179_1_+_112	112	+''',

'foi_duplicate_test': '''chr1\t100\t200
chr1\t100\t200
chr1\t100\t200
''',

'diff_chrom_2': '''chr1	100	400	foi1	.
chr1	70	700	foi2	.
chr2	100	400	foi3	.
chr2	50	800	foi4	.
chr2	1000	1100	foi5	.
chr1	70	800	foi6	.
chr1	1000	1100	foi7	.''',

'same_chrom_2': '''chr1	100	400	foi1	.
chr1	70	700	foi2	.
chr1	100	400	foi3	.
chr1	50	800	foi4	.
chr1	70	800	foi5	.
chr1	1000	1100	foi5	.''',

'gf_chrom12_2': '''chr1	0	100	gf1	.
chr2	0	100	gf1	.''',

'gf_chrom2_2': '''chr2	0	100	gf1	.
''',

'gf_chrom1_2': '''chr1	0	1000	gf1	.
''',

'gf_strand_3': '''chr1	100	200	+
chr1	100	200	.
chr1	100	200	-
''',
'foi_strand_3': '''chr1	100	200	.	.	+
chr1	100	200	.	.	.
chr1	100	200	.	.	-
chr1	300	400	.	.	+
chr1	300	400	.	.	.
chr1	300	400	.	.	-
'''
}

def get_test_data():
	data = {}
	for i,k in feature_data.items():	
		outpath = os.path.join("test",i + ".test.bed")
		with open(outpath, "w") as w:
			w.write(textwrap.dedent(k))
			data[i] = outpath
	return data		

if __name__ == "__main__":
	write_test_data()