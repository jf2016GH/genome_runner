#!/bin/bash
# Modeled after https://superuser.com/questions/302842/resume-rsync-over-ssh-after-broken-connection
#
# Synchronizes UCSC database  and restarts rsync after connection failure
# Usage:
# ./autorsync.sh [dir]
# Where [dir] is the folder to download data

while [ 1 ]
do
    rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/* $1
    if [ "$?" = "0" ] ; then
        echo "rsync completed normally"
        exit
    else
        echo "Rsync failure. Backing off and retrying..."
        sleep 180
    fi
done



