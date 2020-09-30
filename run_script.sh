#!/bin/bash

infile=/Users/rbator01/Box/rebecca_documents/sdmmej_data/PolyGSeq.csv
indir=/Users/rbator01/Box/rebecca_documents/sdmmej_data/
search_radius=30
breakpoint=161
debug=0

reclass=${infile%.csv}_reclassified.csv
deletion=${infile%.csv}_deletion.txt
insertion=${infile%.csv}_insertion.txt

Rscript process_hifibr.R $infile $indir $search_radius $breakpoint $debug

cd deletion/
python SDMMEJDeletionProgram_cli.py -hi $reclass -del $deletion -n $breakpoint -out $indir

cd ../insertion/
Rscript INSERTION_PROGRAM.R $reclass $insertion $indir $breakpoint $search_radius


