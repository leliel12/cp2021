#!/bin/bash
rm -f results/results-headless-peano-*-sorted.dat
for f in results/results-headless-peano-????.dat
do 
	cat ${f} | sort -s -n -k 2,2 | sort -s -n -k 1,1 > ${f%.*}-sorted.dat 
	echo  ${f%.*}-sorted.dat 
done
