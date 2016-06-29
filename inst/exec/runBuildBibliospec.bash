#!/bin/bash

main (){
	find /usr/local/mascot/data -type f \
	| grep "235718\|235705\|235703\|235701\|234943" \
	| grep dat$ | while read i; 
	do 
		~wolski/bin/BiblioSpec/BlibBuild $i `basename $i .dat`.blib; 
	done  
}

main
