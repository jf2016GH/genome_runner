Command line usage of GenomeRunner: SNPs vs. region analysis
===============================================================

Prerequisites: GenomeRunner Web should be installed, and the databases generated. See https://github.com/mdozmorov/genome_runner for more details.

This analysis tests whether a set of SNPs from a region is enriched in the genome annotation data, as compared with randomly selected SNPs from this region.

1. Create a folder to store data necessary for the analysis. Example: ``data.gfs/13qq14_11/``.

2. Create a .bed file with genomic coordinates of the SNPs of interest from the region of interest. Example: "13qq14_11.bed`` (multiple lines).

3. Create a .bed file with genomic coordinates of the region containing all the SNPs of interest, +/- 100bp flanking regions. This region file will be used to extract all the SNPs from the region, to be used as a background. Example: ``13qq14_11_region.bed (single line).

4. Edit ``Makefile.*`` files by changing the name (prefix) of the analysis (e.g., INP = 13qq14_11). Ensure the paths for the background (BKG) and the lists of genomic features (GR) are correct.

5. Run the analysis

.. code-block:: bash

	make -f Makefile.gfs

``data.gfs/13qq14_11/`` - Example of enrichment analysis of selected SNPs using all SNPs from a region as a background

``gfs.hg19.*.txt`` - lists (categories) of genomic features used for the enrichment analysis. See their detailed description in 'Help' section at http://www.genomerunner.org.

``enc.hg19.*.txt`` - lists (categories) of the ENCODE project genomic features. Can be generated using

.. code-block:: bash

	for file in `find /home/mikhail/test_db/grsnp_db/hg19/ENCODE/ -maxdepth 1 ! -path /home/mikhail/test_db/grsnp_db/hg19/ENCODE/ -type d`; do GR=`basename $file`; find $file -type f -name '*.bed.gz' > enc.hg19.$GR.txt; done

``Analysis.Rmd`` - R code for the interpretation of the results. Click ``Analysis.md`` to see an example of the output.