#!/bin/bash
SCRIPT_DIR=$PWD 
INPUT=$1
mkdir -p output
cd $INPUT
for i in $(ls -d */); do
	echo Running parser on: ${i}
	python3 $SCRIPT_DIR/KIRVCFparser.py -i ${i%%/} -o $SCRIPT_DIR/output/${i%%/}_out.csv -j 48 -f .vcf
done