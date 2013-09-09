============
GRSNP Readme
============

------------------------
Installing prerequisites
------------------------

These installation instructions were tested on freshly installed Ubuntu Linux 13.10 (Saucy Salamander). 



.. code-block:: bash

    git clone git@bitbucket.org:mdozmorov/grsnp.git # Clone GenomeRunner SNP repository 
    sudo apt-get install python-setuptools
    sudo apt-get install python-pip
    sudo apt-get install python-dev    
    sudo apt-get upgrade gcc    
    sudo pip install -U cython    
    sudo apt-get install python-numpy python-scipy # Can not be installed via pip install -r     
    # Install R 3.0 (See http://stackoverflow.com/questions/16093331/how-to-install-r-version-3-0)
    # Within R, install packages Hmisc and RColorBrewer by running:
    #    source("http://bioconductor.org/biocLite.R")
    #    biocLite(c("Hmisc", "RCColorBrewer"))
    cd grsnp # Go to grsnp folder, to install the rest of the packages from the file
    sudo pip install -r requirements.txt    
    sudo easy_install simplejson
    # Install GRTK per instructions on https://bitbucket.org/wrenlab/grtk
    python setup.py develop install --user

ToDo: Simplify installation steps
    
There are three important modules in the grsnp package:

1) server # main module
2) dbcreator # A module for creating organism-specific database
3) optimizer # A module to pre-calculate background overlaps with all genomic features


Installing the database
-----------------------

Making the database should be the first step before running GenomeRunner SNP server.

Note that throughout this guide [dir] represents the **full** path to the GRSNP database (i.e [dir] = /home/username/Documents/gr_data)

.. code-block:: bash

    mkdir [dir] # Important: DBcreator should be executed within [dir] folder 
    python -m grsnp.dbcreator -g mm9 -d [dir]

The -d parameter designates where the database is to be installed. This exact path should be passed to grsnp.server.

The -g parameter indicates organism and genome assembly version. Organism-specific genome annotation data are placed in appropriate subfolders, and automatically handled by the grsnp.server.

The DBcreator performs several steps:

* Downloads required files from USCS and place in [dir]/downloads/
* Converts files into standard bed format, sorts them, compresses with *bgzip*.  Unsupported data formats are skipped.
* Places converted files are placed in [dir]/grsnp_db/[organism]/[group]/[tier]/

After creating the database, it is necessary to index bgzipped files with tabix

.. code-block:: bash

    cd [dir]
    find . -type f -name '*.bed.gz' -exec tabix {} \;

ToDo: Automate this into dbcreator.

FAQ
---

* How do I install multiple organism?
  
   * Simply re-run the DBcreator and designate a different organism with the -g parameter.

* Can I run the DBcreator on an existing database?
  
   * Yes, the DBcreator detects which GFs have already been installed.
   
* The DBcreator is taking a long time to run.  Can I 'kill' it?
  
   * Yes, the DBcreator flags partially completed GFs with a '.tmp' extension.  These GFs are not visible to the server.

* Can I download individual GFs?
  
   * Individual GFs can be installed by giving the name to the -f parameter (i.e '-f knownGene' ).

* Is Rsync supported?
  
   * Rsync can be used to mirror the USCS data files .. code-block:: bash
   
       INSERT CODE
   * Before downloading any files, the DBcreator checks the [dir]/downloads/ directory for the required files.
   
Starting and using GenomeRunner SNP
-----------------------------------

GenomeRunner SNP can be started from any folder,but the full path to the database should be provided

.. code-block:: bash
    
    python -m grsnp.server -d [dir]

By default, GenomeRunner SNP uses human hg19 genome annotation data, and 

Pre-calculating GFs and default backgrounds overlap statistics
===============================================================
   
To greatly shorten the enrichment analysis time, grsnp.optimizer should be run on the database.
The number of background SNPs that overlap the GFs (called 'bgs_obs') are obtained once by grsnp.optimizer and stored in a file located at [dir]/grsnp_db/[organism]/bkg_overlaps.gr.
Before running the optimizer, be sure to place some default backgrounds in [dir]/gr_data/custom_data/backgrounds/[organism]/

Continuing from our example above, we can run the following command:

     .. code-block:: bash
     
          python -m grsnp.optimizer -g hg19 -d ~/Documents/gr_data


FAQ
---

* Does the optimizer do all of the organism at once?
  
   * No, the optimizer must be run for each organism
     
* I started the optimizer on an database that already contains the [dir]/grsnp_db/[organism]/bkg_overlaps.gr file.  Can I terminate?
 
   * Yes, you can safely terminate the process.  The file is not modified till the very end of the optimization process.
    
* How is the bkg_overlaps.gr file structured?
  
   [Absolute path to GF file]\t[Absolute path to default background_one]:[bgs_obs]:[n_bgs],[Absolute path to default background_two]:[bgs_obs]:[n_bgs]


Starting the server
===================


A local GRSNP server can be started with the following command.  The -d parameter must point to the directory containing the database created by the dbcreator:

.. code-block:: bash

    python -m grsnp.server -p 8000 -d ~/Documents/gr_data

The server can be access via the following address: 
.. code-block::

    localhost:8081/gr/

The grsnp.server creates several directories in the [/dir] if they do not yet exist.

* [dir]/run_files 
	
   * Stores files that are uploaded to the server by the user and also stores the results of the analysis

The following directories are generated by the sever. Bed files or .gz files should be added to the appropriate directories:

* [dir]/custom_data/gfs/[organism id]
 
   * Stores default genomic features. Sets of default genomic features should be placed in separate folders in this directory.
	
* [dir]/custom_data/fois/[organism id]

   * Stores default feature of interests. Sets of default features of interest files should be placed in separate folders in this directory.

* [dir]/custom_data/backgrounds/[organism id]

   * Stores default backgrounds
 
