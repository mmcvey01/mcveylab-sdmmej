## For renaming the input file
library(tidyverse)

name_file = read.csv("/Users/rbator01/Library/CloudStorage/Box-Box/bioinformatics_research_technology/rt_bioinformatics_consultations/mcvey_lab_rt_bioinformatics/terrence_dna_repair/TestData/test_rename/Iw7R_all_constructs.csv")

name_file$filename1 = paste0("GGAGTACGAA", name_file$X5.UMI,"_", name_file$X3.UMI,"GAACGTTAAC_Final.sam")
name_file$filename2 = paste0( name_file$X5.UMI,"_", name_file$X3.UMI,"_Final.sam")

df = data.frame(name = paste0(name_file$NAME, ".txt"), file1 = name_file$filename1, file2 = name_file$filename2)

nfiles=length(df$name)

is.symlink <- function(paths) isTRUE(nzchar(Sys.readlink(paths), keepNA=TRUE))

path="/Users/rbator01/Library/CloudStorage/Box-Box/bioinformatics_research_technology/rt_bioinformatics_consultations/mcvey_lab_rt_bioinformatics/terrence_dna_repair/HiSeq_CRISPR_Data/all_files/"
newpath="/Users/rbator01/Library/CloudStorage/Box-Box/bioinformatics_research_technology/rt_bioinformatics_consultations/mcvey_lab_rt_bioinformatics/terrence_dna_repair/HiSeq_CRISPR_Data/renamed_files/"
no_events_path="/Users/rbator01/Library/CloudStorage/Box-Box/bioinformatics_research_technology/rt_bioinformatics_consultations/mcvey_lab_rt_bioinformatics/terrence_dna_repair/HiSeq_CRISPR_Data/no_events_files/"

header=c("UMI","CIGAR_STRING","READ_LENGTH","SPLIT_CIGAR_STRING",
         "MATCH_LEFT","MATCH_RIGHT","DISTANCE_FROM_BREAK_LEFT",
         "DISTANCE_FROM_BREAK_RIGHT","DELETION_FROM_LEFT",
         "DELETION_FROM_RIGHT","TOTAL_DELETION","INSERTION_START","INSERTION_END","INSERTION_LENGTH",
         "INSERTED_SEQ","CLASS","ALIGNED_SEQ","READS","MICROHOMOLOGY",
         "MH_Length","NUMBER_OF_ALIGNMENTS","MISMATCH_PERCENTAGE_TO_RECONSTRUCTED")

for ( i in 1:nfiles){
  file1=paste0(path, name_file$filename1[i])
  file2=paste0(path, name_file$filename2[i])
  newfile=paste0(newpath, name_file$NAME[i], ".csv")
  no_events_file=paste0(no_events_path, name_file$NAME[i], ".csv")
  
  if(is.symlink(file1)){
    hifibr_input = read.csv(file1, sep="\t", col.names=header, header=FALSE)
    
    # check that it has non-exact events
    events_for_analysis = hifibr_input %>% dplyr::filter(READS>=10 & CLASS != "exact")
    if (nrow(events_for_analysis) >0){
      write.csv(hifibr_input, newfile, row.names = F, quote=F)
    }else{
      write.csv(hifibr_input, no_events_file, row.names = F, quote=F)
    }
    
  }else if(is.symlink(file2) ){
    hifibr_input = read.csv(file2, sep="\t", col.names=header, header=FALSE)
    
    # check that it has non-exact events
    events_for_analysis = hifibr_input %>% dplyr::filter(READS>=10 & CLASS != "exact")
    if (nrow(events_for_analysis) >0){
      write.csv(hifibr_input, newfile, row.names = F, quote=F)
    }else{
      write.csv(hifibr_input, no_events_file, row.names = F, quote=F)
    }
  }else{
    print(paste0("NOT found ", df$name[i], " ", name_file$filename2[i]))
  }
}
