#!/bin/bash

ls *.chk > chklist
sed -i -e 's/.chk//g' chklist
number=$(wc -l chklist | cut -d" " -f1)
for ii in $(seq 1 $number)
do
	INFILE="`sed -n ${ii}p ./chklist`"
	newzmat -ichk ${INFILE} -oxyz
done
