#!/bin/bash


# Christian Panse <cp@fgcz.ethz.ch>
# 2017-08-17

set -e
set -o pipefail

FASTA="p1875_db10_20170817.fasta"

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

function input(){
	ssh fgcz-176 "cat /srv/localdata/p1644/NL5/NL5_analysis_idx4to5/output/p1644o3482-4_S4.extendedFrags_uniqNB2FC.txt; \
		cat /srv/localdata/p1644/NL5/NL5_analysis_idx4to5/output/p1644o3482-5_S5.extendedFrags_uniqNB2FC.txt"
}


### MAIN
input \
| awk 'NB[$1]{NB[$1]=NB[$1]","$3}
	!NB[$1]{NB[$1]=$3}
	END{
		for(FC in NB){
			printf(">NB%05d NL5idx4-5 FC:%s\n%s\n", ++id, NB[FC], FC)
		}
	}' \
| sed 's/,/,\ /g' \
| cut60 \
> $FASTA \
&&  curl http://fgcz-ms.uzh.ch//FASTA/fgcz_contaminants_20150123.fasta \
>> $FASTA

exit 0
