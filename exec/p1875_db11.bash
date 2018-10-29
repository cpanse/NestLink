#!/bin/bash

# Christian Panse <cp@fgcz.ethz.ch> 2017-09-22

# script to generate the following FASTA
# http://fgcz-ms.uzh.ch/fasta/p1875_db11_20170922.fasta
# Input: txt files containing two columns; 1st is used as FASTAID 2nd as AA sequences
# Output: FASTA
# Usage: p1875_db11.bash textfile > p1875_db11_20170922.fasta

cat $1 \
  | awk '{print ">"$1"\n"$2}' \
  | fgcz_cut60.bash  
  | tee p1875_db11.fasta
