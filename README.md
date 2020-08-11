# sdmmej

This is a repository to work on the error prone dna repair project with the McVey lab at Tufts University.

The script will be run like this on the command line, with a single argument that is the hifibr output file:

Rscript process_hifibr.R TestData_HiFiBR_Output_mod.csv

You need to have the tidyverse and Biostrings library installed.

It generates four output files in a directory like this:

Directory: TestData_HiFiBR_Output_mod_output

Files:
TestData_HiFiBR_Output_mod_reclassified.csv:  same as input, but adds an “ID” column as well as a column for how the sequence was reclassified
TestData_HiFiBR_Output_mod_complex.txt: sequences that could not be reclassified as ins or del
TestData_HiFiBR_Output_mod_insertion.txt: all ins sequences
TestData_HiFiBR_Output_mod_deletion.txt: all del sequences with dashes

I added a few more sequences to the test data, the last three lines of the file TestData_HiFiBR_Output_mod.csv, which are:
1.	Reference sequence
2.	One sequence with <10 reads
3.	One sequence with both SNP and deletion, which ends in the complex category because it can’t be reclassified

