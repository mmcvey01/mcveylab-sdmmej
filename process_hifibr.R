#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Biostrings))

# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("Usage: Rscript process_hifibr.R input_file.csv outdir", call.=FALSE)
} else if (length(args)==2) {
  input_file = args[1]
  out_dir = args[2]
}

#input_file="~/Documents/git/sdmmej/test_data/TestData_HiFiBR_Output_mod.csv"
hifibr_input = read.csv(input_file)

## add an id col and col for reclassified CLASS
hifibr_input <- tibble::rowid_to_column(hifibr_input, "ID")
hifibr_input$CLASS_final = ifelse(hifibr_input$CLASS %in% c('deletion','insertion','exact') , hifibr_input$CLASS, NA)
hifibr_input$CLASS_final = ifelse(hifibr_input$`X._OF_READS` < 10 & hifibr_input$CLASS != 'exact', 'filtered',hifibr_input$CLASS)

## names for output
in_string  = basename(input_file)
out_string = gsub(".csv","",in_string)

## check to make sure only one ref, otherwise exit
ref = hifibr_input %>% dplyr::filter(., CLASS == 'exact')
ref_seq = ref$ALIGNED_SEQ

n_ref = nrow(ref)
if ( n_ref != 1){
  print(paste0("ERROR: Found ", n_ref, " reference sequences in Hifibr output."))
  quit(save = "no", status = 1, runLast = FALSE)
}


# separate sequences and filter for length
hifibr_input_filter_del = hifibr_input %>% dplyr::filter(., `X._OF_READS` > 10 & CLASS == 'deletion')
hifibr_input_filter_ins = hifibr_input %>% dplyr::filter(., `X._OF_READS` > 10 & CLASS == 'insertion')
hifibr_input_filter_complex = hifibr_input %>% dplyr::filter(., `X._OF_READS` > 10 & CLASS == 'complex')

# process deletions
deletions = hifibr_input_filter_del$ALIGNED_SEQ

aligned_deletions = c()

for (seq in deletions){
  sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)
  align <- pairwiseAlignment(seq, ref_seq, substitutionMatrix = sigma, gapOpening = -2,
                                       gapExtension = -8, scoreOnly = FALSE)
  align_string <- toString(align@pattern)
  aligned_deletions = c(aligned_deletions,align_string)
}

#print("we have these aligned deletions")
#print(aligned_deletions)

# process insertions
insertions = hifibr_input_filter_ins$ALIGNED_SEQ

#print("we have these insertions")
#print(insertions)

# process complex
comp_unknown = c()

for (row in 1:nrow(hifibr_input_filter_complex)){
  seq = hifibr_input_filter_complex[row, 'ALIGNED_SEQ']
  id = hifibr_input_filter_complex[row,'ID']
  #print(id)
  align <- pairwiseAlignment(seq, ref_seq, substitutionMatrix = sigma, gapOpening = -2,
                             gapExtension = -8, scoreOnly = FALSE)
  
  ins_width = length(width(insertion(align))[[1]])
  del_width = length(width(deletion(align))[[1]])
  nmismatch = nmismatch(align)
  
  if (ins_width == 0 & del_width > 0 & nmismatch == 0){
    align_string <- toString(align@pattern)
    aligned_deletions = c(aligned_deletions, align_string)
    hifibr_input[hifibr_input$ID == id,]$CLASS_final = 'deletion'
    
    #print("found del")
    #print("ref")
    #print(toString(align@subject))
    #print("seq")
    #print(toString(align@pattern))
  }else if (ins_width > 0 & del_width == 0 & nmismatch == 0){
    insertions = c(insertions, seq)
    hifibr_input[hifibr_input$ID == id,]$CLASS_final = 'insertion'
    
    # print("found ins")
    # print("ref")
    # print(toString(align@subject))
    # print("seq")
    # print(toString(align@pattern))
  }else if (ins_width == 0 & del_width == 0 & nmismatch == 1){
    insertions = c(insertions, seq)
    hifibr_input[hifibr_input$ID == id,]$CLASS_final = 'insertion'
    
    # print("found single mismatch ins")
    # print("ref")
    # print(toString(align@subject))
    # print("seq")
    # print(toString(align@pattern))
  }else{
    # print("found complex")
    # print("ref")
    # print(toString(align@subject))
    # print("seq")
    # print(toString(align@pattern))
    comp_unknown = c(comp_unknown, seq)
    hifibr_input[hifibr_input$ID == id,]$CLASS_final = 'complex'
  }
}

## write out
dir.create(out_dir,showWarnings = FALSE)

## deletions
write(aligned_deletions, file = paste0(out_dir,"/", out_string,"_deletion.txt"),append = FALSE, sep = "\n")

## insertions
insertions = c("RECONSTRUCTED_SEQ",insertions)
write(insertions, file = paste0(out_dir,"/",out_string,"_insertion.txt"),append = FALSE, sep = "\n")

## unknown complex
write(comp_unknown, file = paste0(out_dir,"/", out_string,"_complex.txt"),append = FALSE, sep = "\n")

## output file
write.csv(hifibr_input, file = paste0(out_dir,"/", out_string,"_reclassified.csv"), row.names = FALSE)
