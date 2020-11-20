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

## 
#input_file="~/Box/rebecca_documents/sdmmej_data/rerun_inconsistent_14Oct20/PolyA1Seq_reclassified_inconsistent.csv"
#out_dir="~/Box/rebecca_documents/sdmmej_data/rerun_inconsistent_14Oct20/"
#search_radius=30
#break_location=161
#debug=1

print(paste0("Input file: " ,input_file))
print(paste0("Output dir: ", out_dir))
print(paste0("Search radius: ", search_radius))
print(paste0("Break location: ", break_location))
print(paste0("Debug: ", debug))

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
print(n_ref)

if ( n_ref != 1){
  print(paste0("ERROR: Found ", n_ref, " reference sequences in Hifibr output."))
  quit(save = "no", status = 1, runLast = FALSE)
}

# separate sequences and filter for length
hifibr_input_filter_del = hifibr_input %>% dplyr::filter(.,READS >= 10 & CLASS == 'deletion')
hifibr_input_filter_ins = hifibr_input %>% dplyr::filter(.,READS >= 10 & CLASS == 'insertion')
hifibr_input_filter_complex = hifibr_input %>% dplyr::filter(.,READS >= 10 & CLASS == 'complex')

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

if(debug == 1){
  print("we have these aligned deletions")
  print(aligned_deletions)
}

# process insertions
insertions = hifibr_input_filter_ins$ALIGNED_SEQ

if(debug == 1){
  print("we have these insertions")
  print(insertions)
}
# process complex
comp_unknown = c()

for (row in 1:nrow(hifibr_input_filter_complex)){
  seq = hifibr_input_filter_complex[row, 'ALIGNED_SEQ']
  id = hifibr_input_filter_complex[row,'ID']
  if(debug == 1){
    print(paste0("id is: ",id))
  }
  sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = FALSE)
  align <- pairwiseAlignment(seq, ref_seq, substitutionMatrix = sigma, gapOpening = -2,
                             gapExtension = -8, scoreOnly = FALSE)
  
  if(debug == 1){
    print("ref")
    print(toString(align@subject))
    print("seq")
    print(toString(align@pattern))
    print("insertion")
    print(insertion(align))
    print("deletion")
    print(deletion(align))
    print("mismatch")
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
  
  if(debug==1){
    print("deln")
    print(del_n)
    print("del_pos")
    print(del_pos)
    print("del_r")
    print(del_r)
  }
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
      }
    }else{
      if(debug == 1){
        print(paste0("found del outside search radius at ", del_r))
      }
    }
  }else if (mis_n == 0 & del_n ==0 & ins_n == 1 ){
    if (ins_r < search_radius){
      insertions = c(insertions, seq)
      hifibr_input[hifibr_input$ID == id,]$CLASS_final = 'insertion'
      
      if(debug == 1){
        print("found ins")
      }
    }else{
      if(debug == 1){
        print("found complex (del outside search radius)")
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
      }
    }else{
      if(debug == 1){
        print("found complex (ins outside search radius)")
      }
      comp_unknown = c(comp_unknown, seq)
      hifibr_input[hifibr_input$ID == id,]$CLASS_final = 'complex'
    }
  }else{
    if(debug == 1){
      print("found complex")
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
comp_unknown = c("RECONSTRUCTED_SEQ",comp_unknown)
write(comp_unknown, file = paste0(out_dir,"/", out_string,"_complex.txt"),append = FALSE, sep = "\n")

## output file
write.csv(hifibr_input, file = paste0(out_dir,"/", out_string,"_reclassified.csv"), row.names = FALSE,quote=F)

