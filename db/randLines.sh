#!/usr/bin/env bash
# Randomly shuffles the lines in input file, and generates files of different sizes
# with random content

for i in {10..20}; do ((d=2**i)); shuf snp138.bed | head -n $d > rndsnp138_$d.bed; done

d=0
for i in {1..10}; do ((d=$d+100)); shuf gwascatalog+.bed | head -n $d > rndgwascatalog_$d.bed;done

