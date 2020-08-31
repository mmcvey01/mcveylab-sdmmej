#!/bin/bash

## takes one argument which is the hifi output

input=$1
nick=$2
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

bn=$( basename ${input%.csv} )
results_dir=${DIR}/$( dirname $input )/${bn}_output
mkdir -p ${results_dir}
echo ${results_dir}

hifi_reclass=${results_dir}/${bn}_reclassified.csv
deletion_out=${results_dir}/${bn}_deletion.txt
insertion_out=${results_dir}/${bn}_insertion.txt

echo "------"
echo "Starting Hifibr processing"
#Rscript process_hifibr.R $input $results_dir
echo "Done Hifiber processin"
echo "------"
echo "Starting deletion consistency script"
#python ${DIR}/deletion/SDMMEJDeletionProgram_cli.py -hi ${hifi_reclass} -del ${deletion_out} -n ${nick} -out $results_dir
echo "Done deletion script"
echo "------"
echo "Starting insertion consistency script"
Rscript insertion/INSERTION_PROGRAM.R ${hifi_reclass} ${insertion_out} ${nick} ${results_dir}
echo "Done insertion script"
echo "------"


