Misc. scripts and files for the dbCreator module
================================================

-   autorsync.sh, autorsync1.sh - downloads whole UCSC database using rsync, restarting
    rsync if it breaks
-   blacklisted.txt - genome annotation features that should be ignored
-   extract\_UCSC.py - extract subsets of genome annotation features
    from composite tables into separate files/subfolders
-   Makefile - Automates the ‘extract\_UCSC.py’ to extract most common
    genomic features
-   moveFileLists.sh - Searches selected file names from a file in a
    data folder, and moves them into target folder
-   ucscFilesAll.txt - file list of the UCSC database, sorted by size
-   ucscFilesLarge.txt - largest genome annotation files, to be
    downloaded like
    `` for file in `cat ucscFilesLarge.txt`; do rsync -avzP $file .; done ``
- Makefile.cistrome, Cistrome.*, Epigenome.*, Motif.* - Makefile and prerequisites to download and create Cistrome database for GenomeRunner
- fois* - lists of SNP sets for specific tasks
