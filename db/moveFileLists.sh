#!/usr/bin/env bash
# Searches selected file names from the file (1st argument) in the data folder (2nd argument),
# and moves them out, into defined folder (3rd argument)
# Example: ./moveFileLists.sh ENCODEspecial.txt /home/mikhail/test_db/grsnp_db/ENCODE /home/mikhail/tmp_db/ 

find $2 -type f | fgrep -f $1 | xargs -I{} cp "{}" $3
