#!/bin/bash

# Distributed without any warranties whatsoever.
#
# Christian Panse <cp@fgcz.ethz.ch>
# June 8th 2007

# $HeadURL: http://fgcz-svn.uzh.ch/repos/fgcz/stable/proteomics/bin/fgcz_cut60.bash $
# $Id: fgcz_cut60.bash 7307 2015-03-17 14:29:33Z cpanse $
# $Date: 2015-03-17 15:29:33 +0100 (Tue, 17 Mar 2015) $


test $# -eq 1 || { echo "'$0' requires a input file name (FASTA file) as argument"; exit 1; }
test -f $1 || { echo "'$1'  does not exist"; exit 1; }


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
