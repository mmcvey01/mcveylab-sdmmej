# sdmmej

This is a repository to work on the error prone dna repair project with the McVey lab at Tufts University.

**Dependencies**

Python 2.7

Python libraries: pandas

R >=3.5

R libraries: tidyverse, Biostrings

**Pipeline Script**

This script takes the path to the HiFibr output file as the single command line argument.
Other default arguments are set within the script.

Example usage on test data: 

`sh run_pipeline.sh test_data/polyA1Seq/PolyA1Seq_testdata.csv`

It generates output files in a directory `PolyA1Seq_testdata_output`

- Outputs from Hifibr processing script

`PolyA1Seq_testdata_reclassified.csv` Same format input, but adds an “ID” column as well as a column for how the sequence was reclassified
`PolyA1Seq_testdata_complex.txt` sequences that could not be reclassified as ins or del
`PolyA1Seq_testdata_insertion.txt` all ins sequences
`PolyA1Seq_testdata_deletion.txt` all del sequences with dashes

- Outputs from Deletion script

`PolyA1Seq_testdata_deletion_consistency_log.txt` 
`PolyA1Seq_testdata_deletion_consistency_table.txt`

- Outputs from Insertoin script

`PolyA1Seq_testdatatestdata_insertion_consistency2.csv`
`PolyA1Seq_testdatatestdata_insertion_consistency_long2.csv`
`PolyA1Seq_testdatatestdata_insertion_alignment2.csv`