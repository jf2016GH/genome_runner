Command line usage of GenomeRunner: SNPs vs. region analysis
===============================================================
Prerequisites: GenomeRunner Web should be installed, and the databases generated. See https://github.com/mdozmorov/genome_runner for more details.

``gfs.hg19.*.txt`` - lists (categories) of genomic features used for the enrichment analysis. See their detailed description in 'Help' section at http://www.genomerunner.org.

``SNPs_vs_region`` - Example of selected SNPs enrichment analysis using all SNPs from a region as a background

``enc.hg19.*.txt`` - lists (categories) of the ENCODE project genomic features. Can be generated using

.. code-block:: bash

	for file in `find /home/mikhail/test_db/grsnp_db/hg19/ENCODE/ -maxdepth 1 ! -path /home/mikhail/test_db/grsnp_db/hg19/ENCODE/ -type d`; do GR=`basename $file`; find $file -type f -name '*.bed.gz' > enc.hg19.$GR.txt; done
