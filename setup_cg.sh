#!/usr/bin/env bash

# Requirements:
# - perl
# - python 2 (including headers; python-dev package in ubuntu)
# - R
# - curl
# - GNU tar
# - GNU autotools
# - git
# - gcc

PREFIX=$HOME/.local
mkdir -p $PREFIX/bin
export PATH=$PREFIX:$PREFIX/lib:$PREFIX/bin:$PATH

# install required packages
sudo apt-get install -y build-essential g++
sudo apt-get install -y gdebi-core
sudo apt-get install -y curl
sudo apt-get install -y parallel
sudo apt-get install -y git
sudo apt-get install -y python2.7
sudo apt-get install -y python2.7-dev
sudo apt-get install -y zlib1g-dev # If bedtools errors with fatal error: zlib.h: No such file
# Required by main GR setup
sudo apt-get install -y libreadline-dev
sudo apt-get install -y libpcre3 libpcre3-dev
sudo apt-get install -y liblzma-dev
sudo apt-get install -y libbz2-dev
sudo apt-get install -y libatlas-base-dev
sudo apt-get install -y gfortran

# install R version 3.2.1 for Ubuntu 14.04
which R || {
    wget https://cran.rstudio.com/bin/linux/ubuntu/trusty/r-base-core_3.2.1-4trusty0_amd64.deb
    sudo gdebi -n r-base-core_3.2.1-4trusty0_amd64.deb 
    wget https://cran.rstudio.com/bin/linux/ubuntu/trusty/r-recommended_3.2.1-4trusty0_all.deb
    sudo gdebi -n r-recommended_3.2.1-4trusty0_all.deb 
    wget https://cran.rstudio.com/bin/linux/ubuntu/trusty/r-doc-html_3.2.1-4trusty0_all.deb
    sudo gdebi -n r-doc-html_3.2.1-4trusty0_all.deb 
    wget https://cran.rstudio.com/bin/linux/ubuntu/trusty/r-base_3.2.1-4trusty0_all.deb
    sudo gdebi -n r-base_3.2.1-4trusty0_all.deb 
}

# Versions of software to be installed
declare -A versions
versions=(
    [parallel]=latest
    [kyotocabinet]=1.2.76
    [redis]=2.8.17
    [bedops]=v2.3.0
)

#cd $(mktemp -d)

##########
# Binaries
##########

which bedtools || {
    git clone https://github.com/arq5x/bedtools.git
    cd bedtools
    make
    cp bin/* $PREFIX/bin
    cd -
}

which tabix || {
    git clone https://github.com/samtools/tabix.git
    cd tabix
    make -j
    cp tabix bgzip $PREFIX/bin
    cd -
}

which kchashmgr || {
    base=http://fallabs.com/kyotocabinet/pkg/
    v=${versions[kyotocabinet]}
    curl $base/kyotocabinet-${v}.tar.gz \
        | tar xvzf -
    cd kyotocabinet-${v}
    ./configure --prefix=$PREFIX
    make
    make install
    cd -
}

which realpath || {
    echo -e '#!/usr/bin/env bash\nreadlink -f "$@"' \
        > $PREFIX/bin/realpath
    chmod +x $PREFIX/bin/realpath
}

which redis-server || {
    base=http://download.redis.io/releases/
    v=${versions[redis]}
    curl $base/redis-${v}.tar.gz \
        | tar xvzf -
    cd redis-${v}
    make
    make PREFIX=$PREFIX install
}

which bedops || {
    v=${versions[bedops]}
    curl -L https://github.com/bedops/bedops/releases/download/$v/bedops_linux_x86_64-${v}.tar.bz2 \
        | tar -xjvf - -C $PREFIX
}

cat <<'EOF' |
bedToBigBed
bedGraphToBigWig
bigWigInfo
bigWigSummary
bigWigToBedGraph
bigBedToBed
rowsToCols
EOF
while read p; do
    base=http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
    which $p || {
        curl -o $PREFIX/bin/$p $base/$p
        chmod +x $PREFIX/bin/$p
    }
done

#######################
# GRTK and GenomeRunner
#######################

# install numpy
python -c "import numpy"
if [ $? -gt 0 ]; then
    wget http://ftp.us.debian.org/debian/pool/main/p/python-numpy/python-numpy_1.8.2-2_amd64.deb
    sudo gdebi -n python-numpy_1.8.2-2_amd64.deb 
fi

# install Cython
python -c "import cython"
if [ $? -gt 0 ]; then
    wget http://ftp.us.debian.org/debian/pool/main/c/cython/cython_0.22.1-1_amd64.deb
    sudo gdebi -n cython_0.22.1-1_amd64.deb 
fi

# GenomeRunner branch
branch=shiny

git clone https://github.com/mdozmorov/genome_runner.git
git checkout $branch
cd genome_runner/grtk
python setup.py install --user
cd ..
cd genome_runner
#Rscript installer.R
python setup.py install --user
