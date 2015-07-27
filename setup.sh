#!/usr/bin/env bash
# Installing prerequisites
sudo apt-get install parallel 
sudo apt-get -y install r-base-core 
sudo apt-get install bedtools 
# sudo apt-get install tabix 
sudo apt-get -y install kyotocabinet-utils 
sudo apt-get install realpath
# Ubuntu-specific installation of Python packages. Can be installed using pip install or easy_install
sudo apt-get -y install python-pip python-dev python-cherrypy3 python-numpy python-scipy python-rpy2 python-simplejson python-mako python-beautifulsoup python-celery python-redis
sudo apt-get -y upgrade gcc 
sudo pip install -U cython
sudo apt-get -y install redis-server
sudo pip install celery
sudo pip install flower # Tool to monitor Celery jobs
sudo pip install -U Celery # Helps to solve issue #12, 'module' object has no attribute 'celeryconfiguration'
# Manual download and installation of required binaries
mkdir downloads
cd downloads
wget -N https://github.com/bedops/bedops/releases/download/v2.3.0/bedops_linux_x86_64-v2.3.0.tar.bz2
sudo tar xjvf bedops_linux_x86_64-v2.3.0.tar.bz2 -C /usr/local/
sudo wget -np -R -A "bedToBigBed" -A "bedGraphToBigWig" -A "bigWig*" -A "bigBed*" -N -e robots=off -r -P /usr/local/bin -nd "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/"
sudo wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/rowsToCols
sudo mv rowsToCols /usr/local/bin
# Install tabix v0.2.5
git clone https://github.com/mdozmorov/tabix.git
cd tabix/
make
sudo mv tabix bgzip /usr/local/bin/
cd ..
sudo chmod a+x /usr/local/bin/*
cd ..
sudo rm -r downloads
# Installing R packages
sudo Rscript installer.R
# Installing Genomic Region Tool Kit
cd grtk
sudo python setup.py install
cd ..
# Finally, installing GenomeRunner itself. Keep uncommented only one type of installation
# Developmental mode. Changes made in github-cloned folder are immediately active
sudo python setup.py install develop -d /usr/local/lib/python2.7/dist-packages/
# Standard mode, default. Changes made in github-cloned folder require reinstallation to be active
#sudo python setup.py install
