#!/bin/bash

# Christian Panse <cp@fgcz.ethz.ch> 2017


cat $1 \
| awk '{print ">"$1"\n"$2}' \
| fgcz_cut60.bash  
> p1875_db11.fasta
