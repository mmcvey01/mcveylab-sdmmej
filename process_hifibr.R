#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Biostrings))

# test if there is at least 5 argument: if not, return an error
if (length(args)<5) {
  stop("Usage: Rscript process_hifibr.R input_file.csv outdir search_radius debug", call.=FALSE)
} else if (length(args)==5) {
  input_file = args[1]
  out_dir = args[2]
  search_radius=as.integer(args[3])
  break_location=as.integer(args[4])
  debug=as.integer(args[5])
}

print(paste0("Input file: " ,input_file))
print(paste0("Output dir: ", out_dir))
print(paste0("Search radius: ", search_radius))
print(paste0("Break location: ", break_location))
print(paste0("Debug: ", debug))

## functions

## 
# input_file="~/Documents/git/sdmmej/test_data_2/PolyA1Seq_questions.csv"
# out_dir="~/Documents/git/sdmmej/test_data_2/"
# search_radius=30
# break_location=161
# debug=1

#input_file="~/Documents/git/sdmmej/test_data/TestData_HiFiBR_Output_mod.csv"


hifibr_input = read.csv(input_file)

## add an id col and col for reclassified CLASS
hifibr_input <- tibble::rowid_to_column(hifibr_input, "ID")
hifibr_input$CLASS_final = ifelse(hifibr_input$CLASS %in% c('deletion','insertion','exact') , hifibr_input$CLASS, NA)
hifibr_input$CLASS_final = ifelse(hifibr_input$READS < 10 & hifibr_input$CLASS != 'exact', 'filtered',hifibr_input$CLASS)

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
hifibr_input_filter_del = hifibr_input %>% dplyr::filter(.,READS > 10 & CLASS == 'deletion')
hifibr_input_filter_ins = hifibr_input %>% dplyr::filter(.,READS > 10 & CLASS == 'insertion')
hifibr_input_filter_complex = hifibr_input %>% dplyr::filter(.,READS > 10 & CLASS == 'complex')

# process deletions
deletions = hifibr_input_filter_del$ALIGNED_SEQ

aligned_deletions = c()

for (seq in deletions){
  sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = FALSE)
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

# ### debug
# seq = hifibr_input_filter_complex[2, 'ALIGNED_SEQ']
# id = hifibr_input_filter_complex[2,'ID']
# print(seq)
# print(id)
# sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = FALSE)
# align <- pairwiseAlignment(seq, ref_seq, substitutionMatrix = sigma, gapOpening = -2,
#                            gapExtension = -8, scoreOnly = FALSE)
# 
# ins_info = insertion(align)[[1]]
# del_info = deletion(align)[[1]]
# mis_info = mismatchTable(align)
# 
# ins_n = length(ins_info)
# ins_pos = start(ins_info)
# ins_r = abs(break_location - ins_pos)
# 
# del_n = length(del_info)
# del_pos = start(del_info)
# del_r = abs(break_location - del_pos)
# 
# mis_n = nmismatch(align)
# mis_pos = mis_info$SubjectStart
# mis_r = abs(break_location - mis_pos)
# 
# 
# print(c(ins_n,ins_pos, ins_r))
# print(c(del_n,del_pos, del_r))
# print(c(mis_n,mis_pos, mis_r))
# 
# ## it seems they are always IRanges objects of length 1
# label = ""
# 
# ## check number of deletions
# deletion=deletion(align)[[1]]
# n_del=length(deletion)
# print(n_del)
# # check
# if (n_del >1){
#   label = "complex"  
# }
# 
# # check position of deletions
# pos_del=start(deletion)[1]
# r_del=abs(break_location-pos_del)
# print(r_del)
# if (r_del > search_radius){
#   label = "complex"  
# }
# 
# # check number of insertions
# insertion=insertion(align)[[1]]
# n_ins=length(insertion)
# print(n_ins)
# if (n_ins >1){
#   label="complex"  
# }
# 
# # check position of insertions
# pos_ins=start(insertion)[1]
# print(pos_ins)
# 
# r_ins=abs(break_location-pos_ins)
# print(r_ins)
# # check
# if (r_ins > search_radius){
#   label="complex" 
# }
# 
# # check position of mismatch, assuming it is only one
# pos_mis=mismatchTable(align)$SubjectStart
# r_mis=abs(break_location-pos_mis)
# 
# if (r_mis > search_radius){
#   label="complex" 
# }
# 
# 
# print(label)

### 
for (row in 1:nrow(hifibr_input_filter_complex)){
  seq = hifibr_input_filter_complex[row, 'ALIGNED_SEQ']
  id = hifibr_input_filter_complex[row,'ID']
  if(debug == 1){
    print(id)
  }
  sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = FALSE)
  align <- pairwiseAlignment(seq, ref_seq, substitutionMatrix = sigma, gapOpening = -2,
                             gapExtension = -8, scoreOnly = FALSE)
  
  if(debug == 1){
    print(insertion(align))
    print(deletion(align))
    print(mismatchTable(align))
  }
  
  ins_info = insertion(align)[[1]]
  del_info = deletion(align)[[1]]
  mis_info = mismatchTable(align)
  
  ins_n = length(ins_info)
  ins_pos = start(ins_info)
  ins_r = abs(break_location - ins_pos)
    
  del_n = length(del_info)
  del_pos = start(del_info)
  del_r = abs(break_location - del_pos)
  
  mis_n = nmismatch(align)
  mis_pos = mis_info$SubjectStart
  mis_r = abs(break_location - mis_pos)
  
  # check if it's a single deletion within the search radius
  if (ins_n == 0 & mis_n == 0 & del_n ==1 ){
    if (del_r < search_radius){
      align_string <- toString(align@pattern)
      aligned_deletions = c(aligned_deletions, align_string)
      hifibr_input[hifibr_input$ID == id,]$CLASS_final = 'deletion'
      
      if(debug == 1){
        print("found del")
        print("ref")
        print(toString(align@subject))
        print("seq")
        print(toString(align@pattern))
      }
    }
  }else if (mis_n == 0 & del_n ==0 & ins_n == 1 ){
    if (ins_r < search_radius){
      insertions = c(insertions, seq)
      hifibr_input[hifibr_input$ID == id,]$CLASS_final = 'insertion'
      
      if(debug == 1){
        print("found ins")
        print("ref")
        print(toString(align@subject))
        print("seq")
        print(toString(align@pattern))
      }
    }else{
      if(debug == 1){
        print("found complex (del outside search radius)")
        print("ref")
        print(toString(align@subject))
        print("seq")
        print(toString(align@pattern))
      }
      comp_unknown = c(comp_unknown, seq)
      hifibr_input[hifibr_input$ID == id,]$CLASS_final = 'complex'
    }
  }else if (del_n ==0 & ins_n == 0 & mis_n == 1){
    if (mis_r < search_radius){
      insertions = c(insertions, seq)
      hifibr_input[hifibr_input$ID == id,]$CLASS_final = 'insertion'
      
      if(debug == 1){
        print("found single mismatch ins")
        print("ref")
        print(toString(align@subject))
        print("seq")
        print(toString(align@pattern))
      }
    }else{
      if(debug == 1){
        print("found complex (ins outside search radius)")
        print("ref")
        print(toString(align@subject))
        print("seq")
        print(toString(align@pattern))
      }
      comp_unknown = c(comp_unknown, seq)
      hifibr_input[hifibr_input$ID == id,]$CLASS_final = 'complex'
    }
  }else{
    if(debug == 1){
      print("found complex")
      print("ref")
      print(toString(align@subject))
      print("seq")
      print(toString(align@pattern))
    }
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

