# sdmmej

This is a repository to work on the error prone dna repair project with the McVey lab at Tufts University.

**Installation**  
Installation via Miniconda3 is recommended:

- Configure on the command line to use GitLab on the command line  
`git config --global user.name "tufts username"` (Tufts username is usually 5 letters followed by 2 numbers)  
`git config --global user.email "first.last@tufts.edu"`  

- Download repository  
`git clone https://gitlab.tufts.edu/rbator01/sdmmej.git`  
You will be prompted for password.
Note that you can also download the repository from the web browser if there are problems configuring command line git access.

- Download and install Miniconda3, the following will work on Mac OS  
`curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh --output Miniconda3-latest-MacOSX-x86_64.sh`  
`bash Miniconda3-latest-MacOSX-x86_64.sh`  
follow prompts to accept license and "yes" to running "conda init" on startup  
`source ~/.bash_profile`  

- Create Conda environment called sdmmej  
`conda config --add channels conda-forge`  
`conda config --add channels r`  
`conda config --add channels bioconda`  
`conda create -n sdmmej r-base=4.0.3 r-dplyr r-stringr bioconductor-biostrings python=2.7 pandas`  


- To activate and deactivate you will use  
`conda activate sdmmej`  
`conda deactivate sdmmej`  

**Pipeline Script**

This script takes the path to the HiFibr output file as the single command line argument.
Other default arguments are set within the script.

Example usage on test data: 

`conda activate sdmmej`  
`cd sdmmej-master`  
`sh run_pipeline.sh test_data/polyA1Seq/PolyA1Seq_testdata.csv`  

This generates output files in a directory `PolyA1Seq_testdata_output` with the following output files  

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
