#!/bin/bash


# Christian Panse <cp@fgcz.ethz.ch>
FASTA="p1875_db7_20160704.fasta"

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



# main
cd data_p1875  \
&& for i in 3;
do 
	filename="p1644o2555_$i.extendedFrags_uniqNB2FC.txt"
	cat $filename \
	| grep -v "FlycodeCount" \
	| awk -vfn=$filename 'NF>3{print ">"$NF" "$1" "fn"\n"$3}
		NF<=3{print ">"$1" "fn"\n"$3}'
done \
| sed 's/,//g' \
| cut60 \
> $FASTA

cat /srv/www/htdocs/FASTA/fgcz_contaminants_20150123.fasta \
>>$FASTA
