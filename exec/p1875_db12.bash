#!/bin/bash


# Christian Panse <cp@fgcz.ethz.ch>
# 2017-08-17

set -e
set -o pipefail

FASTA="p1875_db12_20170921.fasta"

function cut60(){
cat $1 \
  | sed 's/\*$//g' \
  | awk '/^>/{print "\n"$0}
	!/^>/{printf("%s",$0)}
  ' \
  | awk '/^>/{print "\n"$0}
	!/^>/{
		for(i=1;i<=length($0);i++){
			printf("%c", substr($0,i,1));
			if (i % 60 == 0) {
				printf("\n");
			}
		}
	}
  ' \
  | egrep "^[A-Z]|^>"
}


for id in 1 2 3 4 5 6; 
do  
  awk -v id=$id '{print id"\t"$0}' \
    < NL6-`echo $id`_S`echo $id`_L001.extendedFrags_uniqNB2FC.txt; 
done \
  | sed 's/,//g'\
  | awk '$2!~/^NB$/{print}' \
  | awk '{printf(">NB%05d_NL6idx%02d FlycodeCount=%d %s\n%s\n", NR, $1, $3, $2, $NF)}' \
  | cut60 \
  > $FASTA \
  &&  curl http://fgcz-ms.uzh.ch//FASTA/fgcz_contaminants_20150123.fasta \
  >> $FASTA

echo $FASTA
