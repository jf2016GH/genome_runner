============
GenomeRunner web server Readme
============

See main documentation on `https://mdozmorov.github.io/grdocs/index.html <https://mdozmorov.github.io/grdocs/index.html>`_

`GenomeRunner web server <http://www.genomerunner.org/>`_ is a tool for functional interpretation of sets of single nucleotide polymorphisms (SNPs) by considering their co-localization with functional/regulatory genome annotation data (epigenomic elements). It is particularly useful for interpretation of functional roles of rare variants and SNPs in non-protein coding regions. An example of GenomeRunner's results can be found in the analysis of Sjogren's syndrome GWAS (`Nature Genetics <http://www.nature.com/ng/journal/v45/n11/full/ng.2792.html>`_
), where it identified RFX5 transcription factor binding site as the most statistically significantly co-localized with the disease-associated SNPs.

GenomeRunner Web calculates enrichment p-values (Chi-squared test) by evaluating whether a set of SNPs co-localizes with regulatory elements more often that could happen by chance. For three or more sets of SNPs, GenomeRunner Web performs the 'epigenomic similarity' analysis by correlating SNP set-specific epigenomic enrichment profiles. Downloadable results are visualized as interactive heatmaps and tables.

This documentation describes the inner workings of GenomeRunner Web interface and command line version, explains the statistics, conventions and terms, and provides the instructions for local installation of GenomeRunner on Linux-based systems.