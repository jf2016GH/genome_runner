============
GRSNP Readme
============

------------------------
Installing prerequisites
------------------------

These installation instructions were tested on freshly installed Ubuntu Linux 13.10 (Saucy Salamander). 



.. code-block:: bash

    git clone git@bitbucket.org:mdozmorov/grsnp.git # Clone GenomeRunner SNP repository 
    sudo apt-get install python-setuptools python-pip python-dev python-numpy python-scipy
    sudo apt-get upgrade gcc
    sudo pip install -U cython   
    # Install R 3.0 (See http://stackoverflow.com/questions/16093331/how-to-install-r-version-3-0)
    # Within R, install packages Hmisc and RColorBrewer by running:
    #    source("http://bioconductor.org/biocLite.R")
    #    biocLite(c("Hmisc", "RCColorBrewer"))
    cd grsnp # Go to grsnp folder, to install the rest of the packages from the file
    sudo pip install -r requirements.txt    
    sudo easy_install simplejson
    # Install GRTK per instructions on https://bitbucket.org/wrenlab/grtk
    python setup.py develop install --user
    # Copy 'frontend' folder and all subfolders into package installation folder.
    # Example: sudo cp -r genome_runner/grsnp/frontend/* .local/lib/python2.7/site-packages/GenomeRunner_SNP-0.1.0-py2.7.egg/grsnp/frontend/

Fixme: Simplify installation steps. Work on eliminating manually copying 'frontend' folder
    
There are three important steps and modules in the grsnp package:

1) dbcreator # A module for creating organism-specific database
2) optimizer # A module to pre-calculate background overlaps with all genomic features
3) server # main module

**dbcreator** # A module for creating organism-specific database
============================================================

Making the database should be the first step before running GenomeRunner SNP server.

Note that throughout this guide [dir] represents the **full** path to the GRSNP database (i.e [dir] = /home/username/Documents/gr_data). Do not use ~ shortcuts.

.. code-block:: bash

    mkdir [dir] # Important: DBcreator should be executed within [dir] folder 
    python -m grsnp.dbcreator -g mm9 -d [dir]

The -d parameter designates where the database is to be installed. This exact path should be passed to grsnp.server. Required.

The -g parameter indicates organism and genome assembly version. Organism-specific genome annotation data are placed in appropriate subfolders, and automatically handled by the grsnp.server. Required.

The ``dbcreator`` performs several steps:

* Downloads required files from USCS and place in [dir]/downloads/
* Converts files into standard bed format, sorts them, compresses with *bgzip*.  Unsupported data formats are skipped, with warning.
* Places converted files are placed in [dir]/grsnp_db/[organism]/[group]. The ENCODE data (tracks beginning with ``wgEncode``) are placed in a separate ``ENCODE`` folder, strictured as [dir]/grsnp_db/[organism]/ENCODE/[data source/type]/[Tier]/[cell type].

After creating the database, it is necessary to index bgzipped files with tabix

.. code-block:: bash

    cd [dir]
    find . -type f -name '*.bed.gz' -exec tabix {} \;

Fixme: Automate this into dbcreator.

Background and custom data
--------------------------
* Custom backgrounds

GenomeRunner uses a 'background', or universe, of all genomic regions for random sampling. For SNPs, a background may be a set of all currently reported SNPs (useful for the analysis of sets of SNPs from Genome-Wide Association Studies), or a set of SNPs on a microarray chip (useful for studies where SNPs selection is limited, such as ImmunoChip, MetaboChip). Such sets of 'background' SNPs should be placed into [dir]/custom_data/backgrounds/[organism]. They are also used to pre-calculate their overlap with each genomic feature by the ``optimizer``.

By default, it is useful to use at least all currently reported as a background. For *Homo Sapiens*, place snp138.bed in the [dir]/custom_data/backgrounds/hg19/ folder. Note that it is important to ensure the end coordinate is larger that the start coordinate, so the processing steps may look like:

   .. code-block:: bash
   
   awk 'BEGIN {OFS="\t"} { if ( $3 <= $2) { print $1, $2, $2+1, $4, $5, $6 } else { print $0 } }' snp138.bed | sort -k1,1 -k2,2n -k3,3n | uniq > snp138+.bed && bgzip snp138+.bed && tabix snp138+.bed.gz 

The logic is, ensure the end coordinate is larger that the start coordinate, sort/unique, block compless, and tabix index the file. 

* Custom features of interest

Sometimes it may be useful to have sets of features of interest readily accessible for the analyses, such as demo sets, or sets of random features. These are placed in subfolders under [dir]/custom_data/fois/[organism]/. The names of the subfolders setve as the descriptions of the sets of fois.

* Custom genomic features

Some genome annotation tracks contain information about different biologically relevant features, lumped together. An example is ``wgEncodeRegTfbsClusteredV3`` track, containing experimentally detected transcription factor binding sites for 161 different transcription factors. The data for each TF can be extracted in separate files using ``extract_UCSC.py`` (see ``db`` subfolder in the source code folder). These files may be placed in [dir]/custom_data/gfs/hg19/tfbsEncode folder, and the 'tfbsEncode' gfs will be accessible through GenomeRunner's interface.

It is a good idea to remove special characters from file names, if any:
   
   .. code-block:: bash
   
   for FILE in *.bed; do mv -v "$FILE" `echo $FILE | tr ' ' '_' | tr -d '[{}(),\!]' | tr -d "\'" | tr '[A-Z]' '[a-z]' | sed 's/_-_/_/g'`;done

and bgzip- and tabix those files for faster processing

   .. code-block:: bash
  
  for file in `find . -type f -name '*.bed'`; do sort -k1,1 -k2,2n -k3,3n $file | uniq > $file"a" && mv $file"a" $file && bgzip $file && tabix $file".gz";done

FAQ
---

* How do I install databases for multiple organism?
  
   * Simply re-run the ``dbcreator`` and designate a different organism with the -g parameter.

* Can I run the ``dbcreator`` on an existing database?
  
   * Yes, the ``dbcreator`` skips GFs that have already been installed.
   
* The ``dbcreator`` is taking a long time to run.  Can I 'kill' it?
  
   * Yes, and you can restart it later. The ``dbcreator`` flags partially completed GFs with a '.tmp' extension.  These GFs are not visible to the server, and will be installed correctly upon next run.

* Can I download individual GFs?
  
   * Individual GFs can be installed by giving the name to the -f parameter (i.e '-f knownGene' ).

* Can I simply download all UCSC data and let the ``dbcreator`` work with it?
  
   * Rsync can be used to mirror the USCS data files. Simply create [dir]/downloads/ folder and execute .. code-block:: bash

      .. code-block:: bash
   
       rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/ .

   Before downloading any files, the ``dbcreator`` checks it they exist in the [dir]/downloads/ folder.
   
**optimizer** # A module to pre-calculate background overlaps with all genomic features
===================================================================================
   
To greatly shorten the enrichment analysis time, the ``optimizer`` should be run on the database. The ``optimizer`` calculates overlap statistics for each genomic feature, that is, how many background regions overlap a genomic feature. The statistics are calculated for each background set of regions ([dir]/custom_data/backgrounds/[organism] folder). These pre-calculated statistics are stored in a file located at [dir]/grsnp_db/[organism]/bkg_overlaps.gr.
Before running the optimizer, be sure to place some default backgrounds in [dir]/custom_data/backgrounds/[organism]/

Continuing from our example above, we can run the following command:

     .. code-block:: bash
     
          python -m grsnp.optimizer -g hg19 -d [dir]

FAQ
---
* Is it necessary to run ``optimizer``?

  * No. If ``bkg_overlaps.gr`` file was not created by the ``optimizer``, GenomeRunner will calculate overlap statistics on the fly. However, calculating overlaps of the background set vs. genomic features on the fly, instead of reading pre-calculated values from the file, takes significant amount of time. So be patient.


* Does the ``optimizer`` do all of the organism at once?
  
  * No, the ``optimizer`` must be run separately for each organism

     
* I started the ``optimizer``, but it takes too long.  Can I terminate?
 
  * Yes, you can safely terminate the process.  The partially completed file bkg_overlaps.gr.tmp will be re-used and appended, when the ``otpimizer`` is restarted.

    
* How is the bkg_overlaps.gr file structured?
  
  * [Absolute path to GF file]\t[Absolute path to default background_one]:[bgs_obs]:[n_bgs],[Absolute path to default background_two]:[bgs_obs]:[n_bgs]

   where [n_bgs] is the total number of regions in the background file, and [bgs_obs] is the number of regions overlapping a genomic feature.


**server** # main module
=========================

GenomeRunner SNP can be started from any folder,but the **full** path to the database should be provided

.. code-block:: bash
    
    python -m grsnp.server -d [dir]

The server can be access via the following address: 

.. code-block::

    localhost:8000/gr/