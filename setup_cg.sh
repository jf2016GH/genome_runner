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

# install requirements
which curl || { 
    apt-get install -y curl
}
which parallel || { 
    apt-get install -y parallel
}
which git || { 
    apt-get install -y git
}
which python2.7 || { 
    apt-get install -y python2.7
}
which python2.7-dev || { 
    apt-get install -y python2.7-dev
}


# GenomeRunner branch
branch=parallel

# Versions of software to be installed
declare -A versions
versions=(
    [parallel]=latest
    [kyotocabinet]=1.2.76
    [redis]=2.8.17
    [bedops]=v2.3.0
)

PREFIX=${PREFIX-$HOME/.local/}
export PATH=$PREFIX/bin:$PATH
mkdir -p $PREFIX/bin

#cd $(mktemp -d)

##########
# Binaries
##########

which parallel || {
    v=${versions[parallel]}
    curl http://ftp.gnu.org/gnu/parallel/parallel-${v} \
        | tar xjvf -
    cd parallel-${v}
    ./configure --prefix=$PREFIX
    make install
    cd -
}

which bedtools || {
    git clone https://github.com/arq5x/bedtools.git
    cd bedtools
    make -j
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
    make -j
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
    make -j hiredis jemalloc linenoise lua
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

# git clone https://github.com/mdozmorov/genome_runner.git
# git checkout $branch
# cd genome_runner/grtk
# python setup.py install --user
# cd ..
# cd genome_runner
# Rscript installer.R
# python setup.py install --user
