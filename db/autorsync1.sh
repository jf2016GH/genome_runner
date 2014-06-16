#!/bin/bash
# Modeled after http://www.castledragmire.com/Posts/Automatically_resuming_rsync
#
# Synchronizes UCSC database  and restarts rsync after connection failure
# Usage:
# ./autorsync.sh [dir]
# Where [dir] is the folder to download data

export Result=1; #This will hold the result of the rsync. Set to 1 so the first loop check will fail.
while [ $Result -ne 0 ]; do #Loop until rsync result is successful
  echo "STARTING ($Result) @" `date`; #Inform the user of the time an rsync is starting and the last rsync failure code
  rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/* $1
  Result=$?; #Store the result of the rsync
  sleep 1; #This is an optional 1 second timeout between attempts
done
