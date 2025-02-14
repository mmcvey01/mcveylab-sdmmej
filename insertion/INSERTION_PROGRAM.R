# 08/05/2020
# Script to find MOST templated/SD-MMEJ consistent insertions.

# Currently, it will only give you the smallest motif repeat of +1 nucleotide on either side of the insertion. This program only gives you a start with the insertions. You need to manually go through and see if the repeat motif can be longer
# Also, it cannot analyize insertions where the insertions are "ttattat", i.e. fill in synthesis of the ttat overhangs. You will have to manually analyize those.

# Start by defining the sequence of the region to search. Break up the sequence into "left" and "right" of the break. If the break leaves overhangs, include overhangs on both the left and right sequences.
# also define the distance around the break you wish to search for SD-MMEJ consistent repeats.

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressWarnings(suppressMessages({library(tidyverse)
  library(Biostrings)
  library(stringr)
  library(dplyr)
}))

#test if there is at least one argument: if not, return an error
if (length(args)<5) {
  print(length(args))
  stop("Usage: Rscript INSERTION_PROGRAM.R hifibr_reclassified.csv insertion_file.txt nick_location outdir", call.=FALSE)
} else if (length(args)==5) {
  hifi_in = args[1]
  insertion_in = args[2]
  out_dir = args[3]
  nick = as.integer(args[4])
  search_radius=as.integer(args[5])

}

# hifi_in="/Users/rbator01/Library/CloudStorage/Box-Box/bioinformatics_research_technology/rt_bioinformatics_consultations/mcvey_lab_rt_bioinformatics/terrence_dna_repair/HiSeq_CRISPR_Data/renamed_files/R3_output/R3_reclassified.csv"
# insertion_in="/Users/rbator01/Library/CloudStorage/Box-Box/bioinformatics_research_technology/rt_bioinformatics_consultations/mcvey_lab_rt_bioinformatics/terrence_dna_repair/HiSeq_CRISPR_Data/renamed_files/R3_output/R3_insertion.txt"
# outdir="/Users/rbator01/Library/CloudStorage/Box-Box/bioinformatics_research_technology/rt_bioinformatics_consultations/mcvey_lab_rt_bioinformatics/terrence_dna_repair/HiSeq_CRISPR_Data/renamed_files/R3_output/"
# nick=161
# search_radius=30

a<-read.csv(insertion_in)

if(nrow(a) == 0){
  stop("\r No insertions in file")
}
## get reference
hifibr_input = read.csv(hifi_in,header=T)


ref = hifibr_input %>% dplyr::filter(., CLASS_final == 'exact')
n_ref = nrow(ref)

if ( n_ref != 1){
  print(paste0("ERROR: Found ", n_ref, " reference sequences in Hifibr output."))
  quit(save = "no", status = 1, runLast = FALSE)
}


ref = ref$ALIGNED_SEQ

n <- search_radius  # number of bases to the left and right of the break you want to search of repeated motifs
p <- 10 # number of bases to the left and right of the initial repeat motif to search for homology

### for naming the output files
plasmid <- gsub("_reclassified.csv","",basename(hifi_in))
type <- gsub(plasmid, "", basename(insertion_in))
type <- gsub(".txt", "", type)
type <- gsub("_", "", type)

L = substr(ref,0,nick)
print("left of nick")
print(L)

l=length(ref)
R <- substr(ref,nick+1,nchar(ref)) 
print("right of nick")
print(R)

l <- nchar(L) # number of nucleotides of the left hand sequence
r <- nchar(R) # number of nucleotides of the right hand sequence

k1 <- 30 # how far you want to cut back to search, this needs to to be adjusted based on how large the deletions are, but i think this covers up to 30 bp of deletion on either side.
k2 <- 30

sL1 <- substring(L, 1, (l-k1-1):l)
sR1 <- substring(R, 1:(r-k2-1),r)

a2=NULL    # create empty vector to insert left del boundary
for (i in a[, 1]){
  lb <- str_locate(as.character(i), sL1)
  lb <- na.omit(lb)
  lbb <- lb[nrow(lb),2]
  if (length(lbb) == 0){
    print(paste0("No match found left of boundary for sequence ", i))
    a2[i] = 0
  }else{
    a2[i] = lbb
  }
}
a3=NULL    # create empty vector to insert RIGHT del boundary
for (i in a[, 1]){
  rb <- str_locate(as.character(i), sR1)
  rb <- na.omit(rb)
  rbb <- rb[1,1]
  if (length(rbb) == 0){
    print(paste0("No match found reft of boundary for sequenc ", i))
    a3[i] = 0
  }else{
    a3[i] = rbb
  }
}


a4 <- cbind(as.data.frame(a2), as.data.frame(a3))    # combine left and right del boundary
a5 <- cbind(a, a4)    # combine seq and del boundary
names(a5) <- c("RECONSTRUCTED_SEQ","left_del", "right_del")    # rename columns for ease

del_seq <- paste(substring(a5$RECONSTRUCTED_SEQ, first = 1, last = a5$left_del),    # create sequence without the insertion sequence
                 substring(a5$RECONSTRUCTED_SEQ, first = a5$right_del, last = nchar(as.character(a5$RECONSTRUCTED_SEQ))),sep = "") 

a5$del_seq <- del_seq
ins <- substring(a5$RECONSTRUCTED_SEQ, first = a5$left_del+1, last = a5$right_del-1)    # extract inserted sequence
a5$insertion <- ins    # Now I have insertion added to table
a5$plasmid <- plasmid
a5$ID <- paste(a5$plasmid, 1:nrow(a5), sep="-")

# SOMETHING IS WRONG WITH NO.BUTT!!!! 06282019
no.butt <- a5
#[-which(a5$insertion == ""), ] # dont want to have the insertion sequences that are all fucky, i.e. have ttat 
# overhang insert define the region around the break site you want to look for homology. NOTE: value is one less than
# desired...dont worry

a6 <- data.frame(ID=numeric() ,x=numeric() , y=numeric() , stringsAsFactors = FALSE)
x <- 10
for (i in 1:nrow(no.butt)){
  ins <- substring(no.butt[i,1], first = no.butt[i,2], last = no.butt[i,3])    # create insertion sequence region to search
  dr1 <- as.data.frame(str_locate(no.butt[i,4], ins))    # search deletion sequence for string
  if(!is.na(dr1[1,1])){
    dr2 <- as.data.frame(str_locate_all(no.butt[i,4], ins))
    a7 <- cbind(no.butt[i,7], dr2)
    a6 <- rbind(a6,a7)
    while (nrow(dr2)!=0){
      for (j in 1:p) {
        ins <- substring(no.butt[i,1], first = no.butt[i,2]-j, last = no.butt[i,3])
        dr3 <- as.data.frame(str_locate(no.butt[i,4], ins))
        if (!is.na(dr3[1,1])){
          dr2 <- as.data.frame(str_locate_all(no.butt[i,4], ins))
          a7 <- cbind(no.butt[i,7], dr2)
          a6 <- rbind(a6,a7)
        } else {
          for (k in 1:p){
            ins <- substring(no.butt[i,1], first = no.butt[i,2]-j+1, last = no.butt[i,3]+k)
            dr3 <- as.data.frame(str_locate(no.butt[i,4], ins))
            if (!is.na(dr3[1,1])){
              dr2 <- as.data.frame(str_locate_all(no.butt[i,4], ins))
              a7 <- cbind(no.butt[i,7], dr2)
              a6 <- rbind(a6,a7)
            }
            else{
              dr2 <- data.frame()
              break
            }
          }
          dr2 <- data.frame()
          break
        }
      }
      if(j==x){
        break
      }
      if(k==x){
        break
      }
    }
  }else {
    a7 <- cbind(no.butt[i,7], dr1)
    a6 <- rbind(a6,a7)
  }
}
# Now search on the right side first
a8 <- data.frame(ID=numeric() ,x=numeric() , y=numeric() , stringsAsFactors = FALSE)
x <- 10   # define the distance to the left and right you wish to search
for (i in 1:nrow(no.butt)){
  ins <- substring(no.butt[i,1], first = no.butt[i,2], last = no.butt[i,3])    # create insertion sequence region to search
  dr1 <- as.data.frame(str_locate(no.butt[i,4], ins))    # search deletion sequence for string
  if(!is.na(dr1[1,1])){
    dr2 <- as.data.frame(str_locate_all(no.butt[i,4], ins))
    a7 <- cbind(no.butt[i,7], dr2)
    a8 <- rbind(a8,a7)
    while (nrow(dr2)!=0){
      for (j in 1:p) {
        ins <- substring(no.butt[i,1], first = no.butt[i,2], last = no.butt[i,3]+j)
        dr3 <- as.data.frame(str_locate(no.butt[i,4], ins))
        if (!is.na(dr3[1,1])){
          dr2 <- as.data.frame(str_locate_all(no.butt[i,4], ins))
          a7 <- cbind(no.butt[i,7], dr2)
          a8 <- rbind(a8,a7)
        } else {
          for (k in 1:p){
            ins <- substring(no.butt[i,1], first = no.butt[i,2]-k, last = no.butt[i,3]+j-1)
            dr3 <- as.data.frame(str_locate(no.butt[i,4], ins))
            if (!is.na(dr3[1,1])){
              dr2 <- as.data.frame(str_locate_all(no.butt[i,4], ins))
              a7 <- cbind(no.butt[i,7], dr2)
              a8 <- rbind(a8,a7)
            }
            else{
              dr2 <- data.frame()
              break
            }
          }
          dr2 <- data.frame()
          break
        }
      }
      if(j==x){
        break
      }
      if(k==x){
        break
      }
    }
  }else {
    a7 <- cbind(no.butt[i,7], dr1)
    a8 <- rbind(a8,a7)
  }
}

revC <- data.frame(ID=numeric() ,x=numeric() , y=numeric() , stringsAsFactors = FALSE)
x <- 10   # define the distance to the left and right you wish to search
for (i in 1:nrow(no.butt)){
  ins2 <- substring(no.butt[i,1], first = no.butt[i,2], last = no.butt[i,3])    # create insertion sequence region to search
  dnains2 <- lapply(ins2, DNAString)
  revins2 <- lapply(dnains2, reverseComplement)
  revins2 <- lapply(revins2, toString)
  dr1 <- as.data.frame(str_locate(no.butt[i,4], as.character(revins2)))    # search deletion sequence for string
  if(!is.na(dr1[1,1])){
    dr2 <- as.data.frame(str_locate_all(no.butt[i,4], as.character(revins2)))
    a7 <- cbind(no.butt[i,7], dr2)
    revC <- rbind(revC,a7)
    while (nrow(dr2)!=0){
      for (j in 1:p) {
        ins2 <- substring(no.butt[i,1], first = no.butt[i,2]-j, last = no.butt[i,3])
        dnains2 <- lapply(ins2, DNAString)
        revins2 <- lapply(dnains2, reverseComplement)
        revins2 <- lapply(revins2, toString)
        dr3 <- as.data.frame(str_locate(no.butt[i,4], as.character(revins2)))
        if (!is.na(dr3[1,1])){
          dr2 <- as.data.frame(str_locate_all(no.butt[i,4], as.character(revins2)))
          a7 <- cbind(no.butt[i,7], dr2)
          revC <- rbind(revC,a7)
        } else {
          for (k in 1:p){
            ins2 <- substring(no.butt[i,1], first = no.butt[i,2]-j+1, last = no.butt[i,3]+k)
            dnains2 <- lapply(ins2, DNAString)
            revins2 <- lapply(dnains2, reverseComplement)
            revins2 <- lapply(revins2, toString)
            dr3 <- as.data.frame(str_locate(no.butt[i,4], as.character(revins2)))
            if (!is.na(dr3[1,1])){
              dr2 <- as.data.frame(str_locate_all(no.butt[i,4], as.character(revins2)))
              a7 <- cbind(no.butt[i,7], dr2)
              revC <- rbind(revC,a7)
            }
            else{
              dr2 <- data.frame()
              break
            }
          }
          dr2 <- data.frame()
          break
        }
      }
      if(j==x){
        break
      }
      if(k==x){
        break
      }
    }
  }else {
    a7 <- cbind(no.butt[i,7], dr1)
    revC <- rbind(revC,a7)
  }
}
# search from right to left this time
revC2 <- data.frame(ID=numeric() ,x=numeric() , y=numeric() , stringsAsFactors = FALSE)
x <- 10   # define the distance to the left and right you wish to search
for (i in 1:nrow(no.butt)){
  ins2 <- substring(no.butt[i,1], first = no.butt[i,2], last = no.butt[i,3])    # create insertion sequence region to search
  dnains2 <- lapply(ins2, DNAString)
  revins2 <- lapply(dnains2, reverseComplement)
  revins2 <- lapply(revins2, toString)
  dr1 <- as.data.frame(str_locate(no.butt[i,4], as.character(revins2)))    # search deletion sequence for string
  if(!is.na(dr1[1,1])){
    dr2 <- as.data.frame(str_locate_all(no.butt[i,4], as.character(revins2)))
    a7 <- cbind(no.butt[i,7], dr2)
    revC2 <- rbind(revC2,a7)
    while (nrow(dr2)!=0){
      for (j in 1:p) {
        ins2 <- substring(no.butt[i,1], first = no.butt[i,2]-j, last = no.butt[i,3])
        dnains2 <- lapply(ins2, DNAString)
        revins2 <- lapply(dnains2, reverseComplement)
        revins2 <- lapply(revins2, toString)
        dr3 <- as.data.frame(str_locate(no.butt[i,4], as.character(revins2)))
        if (!is.na(dr3[1,1])){
          dr2 <- as.data.frame(str_locate_all(no.butt[i,4], as.character(revins2)))
          a7 <- cbind(no.butt[i,7], dr2)
          revC2 <- rbind(revC2,a7)
        } else {
          for (k in 1:p){
            ins2 <- substring(no.butt[i,1], first = no.butt[i,2]-j+1, last = no.butt[i,3]+k)
            dnains2 <- lapply(ins2, DNAString)
            revins2 <- lapply(dnains2, reverseComplement)
            revins2 <- lapply(revins2, toString)
            dr3 <- as.data.frame(str_locate(no.butt[i,4], as.character(revins2)))
            if (!is.na(dr3[1,1])){
              dr2 <- as.data.frame(str_locate_all(no.butt[i,4], as.character(revins2)))
              a7 <- cbind(no.butt[i,7], dr2)
              revC2 <- rbind(revC2,a7)
            }
            else{
              dr2 <- data.frame()
              break
            }
          }
          dr2 <- data.frame()
          break
        }
      }
      if(j==x){
        break
      }
      if(k==x){
        break
      }
    }
  }else {
    a7 <- cbind(no.butt[i,7], dr1)
    revC2 <- rbind(revC2,a7)
  }
}
a10 <- rbind(a6, a8)
a10$motif_length <- a10[,3]-a10[,2]
a10 <- a10[order(a10$motif_length, decreasing=TRUE),]
a10$Lunicorn <- paste(a10[,1], a10[,2], sep="-") 
a10$Runicorn <- paste(a10[,1], a10[,3], sep="-") 
a11 <- a10[!duplicated(a10[,5]),]
a12 <- a11[!duplicated(a11[,6]),]
revCm <- rbind(revC, revC2)
revCm$motif_length <- revCm[,3]-revCm[,2]
revCm <- revCm[order(revCm$motif_length, decreasing=TRUE),]
revCm$Lunicorn <- paste(revCm[,1], revCm[,2], sep="-") 
revCm$Runicorn <- paste(revCm[,1], revCm[,3], sep="-") 
revCm2 <- revCm[!duplicated(revCm[,5]),]
revCm3 <- revCm2[!duplicated(revCm2[,6]),]
a6 <- a12[,1:3]
revC <- revCm3[,1:3]
names(a6) <- c("ID", "DR_START", "DR_END")
names(revC) <- c("ID", "RC_START", "RC_END")
a8 <- merge(a6, no.butt, by="ID")
a8$left_yes <- ifelse(a8$DR_START>(a8$left_del-n), "1", "0")
a8$right_yes <- ifelse(a8$DR_END<(a8$right_del+n), "1", "0")
a8$sum <- (as.numeric(a8$left_yes)+as.numeric(a8$right_yes))
a9 <- subset(a8, a8$sum==2)
a9 <- rbind(a9, subset(a8, is.na(a8$sum)))
a6 <- a9[,1:3]
revC2 <- merge(revC, no.butt, by="ID")
revC2$left_yes <- ifelse(revC2$RC_END>(revC2$left_del-n), "1", "0")
revC2$right_yes <- ifelse(revC2$RC_START<(revC2$right_del+n), "1", "0")
revC2$sum <- (as.numeric(revC2$left_yes)+as.numeric(revC2$right_yes))
revC3 <- subset(revC2, revC2$sum==2)
revC3 <- rbind(revC3, subset(revC2, is.na(revC2$sum)))
revC <- revC3[,1:3]
a6$consistency  <- ifelse(is.na(a6$DR_START), "0", "1")    # to add consistency
revC$consistency  <- ifelse(is.na(revC$RC_START), "0", "1")     # to add consistency
master <- merge(a6, revC, by="ID", all=TRUE)    # combine reverse complement with direct repeat
master$consistency.x <- ifelse(is.na(master$DR_START), "0", "1")
master$consistency.y <- ifelse(is.na(master$RC_START), "0", "1")
master$sum <- (as.numeric(master$consistency.x)+as.numeric(master$consistency.y))
master$consistency <- ifelse(master$sum==0, "FALSE", "TRUE")    # consistency for either direct repeat/reverse complement
master <- cbind(master[,1:3], master[,5:6], master[,9])
names(master) <- c("ID", "DR_START", "DR_END","RC_START", "RC_END", "consistency")
master2 <- merge(master, a5, by="ID", all.y=TRUE) # add to original data set!
master2$DRmotif_length  <- master2$DR_END-master2$DR_START+1    # add repeat motif length for direct repeat
master2$RCmotif_length <- master2$RC_END-master2$RC_START+1    # add repeat motif length for reverse complement
master2[,2] <- ifelse(master2$DRmotif_length>=4,master2[,2] , NA)    # consistency for either direct repeat/reverse complement
master2[,4] <- ifelse(master2$RCmotif_length>=4,master2[,4] , NA)
master2[,6] <- ifelse(!is.na(master2[,2]), "TRUE", ifelse(!is.na(master2[,4]), "TRUE", "FALSE" ))
# Code below workds to give you "aligned" templates
for (i in 1:nrow(master2)){
  if(!is.na(master2[i,2])){
    test <- substring(master2[i,10], first = master2[i,2], last = master2[i,3])
    test2 <- ifelse((as.numeric(master2[i,2])>as.numeric(master2[i,8])), (as.numeric(master2[i,2])+as.numeric(nchar(master2[i,11]))-1), (as.numeric(master2[i,2])-1))
    test3 <- rep.int("-", times=test2)
    test4 <- paste(test3, collapse="")
    test5 <- paste(test4, test, sep="")
    test6 <- rep.int("-", times=(as.numeric(nchar(as.character(master2[i,7]))-nchar(test5))))
    test7 <- paste(test6, collapse="")
    test8 <- paste(test5, test7, sep="")
    master2[i,15]=test8    
  }
  else {
    master2[i,15]="0"
  }
}
# For reverse complement
for (i in 1:nrow(master2)){
  if(!is.na(master2[i,4])){
    test <- tolower(substring(master2[i,10], first = master2[i,4], last = master2[i,5]))    # make reverse complement repeats lowercase
    test2 <- ifelse((as.numeric(master2[i,4])>as.numeric(master2[i,8])), (as.numeric(master2[i,4])+as.numeric(nchar(master2[i,11]))-1), (as.numeric(master2[i,4])-1))
    test3 <- rep.int("-", times=test2)
    test4 <- paste(test3, collapse="")
    test5 <- paste(test4, test, sep="")
    test6 <- rep.int("-", times=(as.numeric(nchar(as.character(master2[i,7]))-nchar(test5))))
    test7 <- paste(test6, collapse="")
    test8 <- paste(test5, test7, sep="")
    master2[i,16]=test8    
  }
  else {
    master2[i,16]="0"
  }
} 
colnames(master2)[15:16] <- c("Loop-out", "Snap-back")
test <- reshape(master2, 
                varying=c("RECONSTRUCTED_SEQ", "Loop-out", "Snap-back"),
                idvar="ID",
                v.names="insertion_alignment",
                timevar="mechanism",
                times=c("seq", "Loop-out", "Snap-back"),
                new.row.names=1:10000,
                direction="long")
test <- test[order(test$ID),]
test$unicorn <- paste(test[,1], test[,14], test[,15], sep="-") 
test2 <- test[!duplicated(test[,16]),]
# save it, before making a slightly cleaner version
# Make pretty, two columns, ID, and aligned seq
pretty <- cbind(ID=test2[,1], as.data.frame(test2[,15]), as.data.frame(test2[,14])) 
colnames(pretty)[2:3] <- c("insertion_alignment", "mechanism")

output1 <- paste0(out_dir, "/", plasmid, "_",type,"_insertion_consistency2.csv")
output2 <- paste0(out_dir, "/", plasmid, "_",type,"_insertion_consistency_long2.csv")
output3 <- paste0(out_dir, "/", plasmid, "_",type,"_insertion_alignment2.csv")
write.csv(master2, output1)
write.csv(test2, output2)
write.csv(pretty, output3)

