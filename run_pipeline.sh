#!/bin/bash

## This script takes one input: the path to the HiFibr output file
## Example usage on test data: sh run_pipeline.sh test_data/polyA1Seq/PolyA1Seq_testdata.csv
if [ "$#" -ne 1 ]; then
    echo "Error: please supply the path to the HiFibr input file as a command line argument"
fi


input=$1

## These parameters are the defaults, can be changed in this script if needed
search_radius=30
breakpoint=161
debug=0
##


bn=$( basename ${input%.csv} )
#results_dir=$( pwd )/$( dirname $input )/${bn}_output
results_dir=$( dirname $input )/${bn}_output
mkdir -p ${results_dir} ${results_dir}/table_outputs ${results_dir}/plots

hifi_reclass=${results_dir}/${bn}_reclassified.csv
deletion_out=${results_dir}/${bn}_deletion.txt
insertion_out=${results_dir}/${bn}_insertion.txt
complex_out=${results_dir}/${bn}_complex.txt

echo ${hifi_reclass}

echo "------"
echo "Starting to process HiFibr output file"
echo "------"

# Rscript process_hifibr.R $input $results_dir $search_radius $breakpoint $debug
# 
# echo "------"
# echo "Done Hifiber processing"
# echo "------"
# echo "Starting deletion consistency script, see log ${bn}_deletion.log"
# 
# python deletion/SDMMEJDeletionProgram_cli.py -hi ${hifi_reclass} -del ${deletion_out} -n $breakpoint -out $results_dir > ${results_dir}/${bn}_deletion.log
# 
# echo "------"
# echo "Done deletion script"
# echo "------"
# echo "------"
# echo "Starting insertion consistency script on insertions"
# echo "------"
# 
# Rscript insertion/INSERTION_PROGRAM.R ${hifi_reclass} ${insertion_out} $results_dir $breakpoint $search_radius
# 
# echo "------"
# echo "Done insertion script on insertions"
# echo "------"
# 
# echo "------"
# echo "Starting insertion consistency script on complex"
# echo "------"
# 
# Rscript insertion/INSERTION_PROGRAM.R ${hifi_reclass} ${complex_out} $results_dir $breakpoint $search_radius
# 
# echo "------"
# echo "Done insertion script on complex"
# echo "------"
# 
# ## Run plotting script

Rscript SDMMEJ_Plots.R $results_dir $bn

echo "------"
echo "Done plotting script"
echo "------"

echo "Results will be located in ${results_dir}" 
