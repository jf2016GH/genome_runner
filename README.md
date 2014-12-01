--------------------------------
  GenomeRunner web Readme
  --------------------------------

See the main documentation on
[https://mdozmorov.github.io/grdocs/index.html](https://mdozmorov.github.io/grdocs/index.html)

[GenomeRunner web](http://www.genomerunner.org/) is a tool for  investigation of potential regulatoy impact of sets of single nucleotide polymorphisms (SNPs) by considering their
co-localization with functional/regulatory genome annotation data (regulatory datasets). It is particularly useful for the interpretation of functional roles of sets of rare variants and SNPs in non-protein coding regions. An example of GenomeRunner’s results can be found in the analysis of Sjogren’s syndrome GWAS ([Nature Genetics](http://www.nature.com/ng/journal/v45/n11/full/ng.2792.html) ), where it identified RFX5 transcription factor binding site as the most statistically significantly co-localized with the set of disease-associated SNPs.

GenomeRunner web calculates enrichment p-values (Chi-squared test) by evaluating whether a set of SNPs co-localizes with regulatory datasets more often that could happen by chance. For three or more sets of SNPs, GenomeRunner web performs the ‘regulatory similarity’ analysis by correlating SNP set-specific regulatory enrichment profiles. Downloadable results are visualized as interactive heatmaps and tables.
