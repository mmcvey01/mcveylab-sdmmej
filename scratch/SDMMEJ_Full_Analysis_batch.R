# ========================================================== Initial Analysis for SD-MMEJ Pipeline Data =======================================================================
library(seqinr)
library(ggplot2)
library(grid)
library(gtable)
library(stringr)
library(dplyr)

## input files
# a<-read.csv("Iw7Cas9_HiFiBR_reclassified.csv") # VARIABLE
# in_template_name="Iw7"
# in_consistency_log = "Iw7_python_consistency.csv"
# in_ins_consistency = "Iw7Cas9_HiFiBR_insertion_insertion_consistency2.csv" 
# in_comp_consistency = "Iw7Cas9_HiFiBR_complex_insertion_consistency2.csv"
# in_break_curated="Iw7break_curated.csv"
# 
# ## output files
# out_cc ="Iw7_combined_curated.csv"
# out_dmerge="Iw7_del_merge_table.csv"
# out_all="IW7_all_SD-MMEJ_consistency.csv"
# out_break="Iw7break.csv"
# out_ins_consistency="Iw7_insertion_consistent_analysis.csv"
# out_del_for_template_plot="Iw7_del_data_for_template_plot.csv"
# out_del_compressed_temp_plot_mech="Iw7_del_data_for_compressed_temp_plot-mech.csv"
# out_del_compressed_temp_plot="Iw7_del_data_for_compressed_temp_plot.csv"
# out_ins_template_plot_data="Iw7_ins_template_plot_data.csv"
# out_ins_template_plot_data_mechanism="Iw7_ins_template_plot_data_mechanism.csv"
# out_consistency_true_ins="Iw7_SD-MMEJ_consistency_true_ins_del2.csv"

## input files
setwd('/Users/rbator01/Box/git/sdmmej/test_data/polyA1Seq/PolyA1Seq_testdata_output/')
in_template_name="PolyA1Seq_testdata"

a<-read.csv(paste0(in_template_name, "_reclassified.csv")) # VARIABLE
in_consistency_log = paste0(in_template_name, "_deletion_consistency_log.txt")
in_ins_consistency = paste0(in_template_name, "_insertion_insertion_consistency2.csv" )
in_comp_consistency = paste0(in_template_name, "_complex_insertion_consistency2.csv")
in_break_curated=paste0(in_template_name, "_break.csv") #"Iw7break_curated.csv" QUESTION IF CURATION IS NECESSARY

## output files
outfolder="output/"
out_cc =paste0(outfolder,in_template_name, "_combined_curated.csv")
out_dmerge=paste0(outfolder,in_template_name, "_del_merge_table.csv")
out_all=paste0(outfolder,in_template_name, "_all_SD-MMEJ_consistency.csv")
out_break=paste0(outfolder,in_template_name, "_break.csv")
out_ins_consistency=paste0(outfolder,in_template_name, "_insertion_consistent_analysis.csv")
out_del_for_template_plot=paste0(outfolder,in_template_name, "_del_data_for_template_plot.csv")
out_del_compressed_temp_plot_mech=paste0(outfolder,in_template_name, "_del_data_for_compressed_temp_plot-mech.csv")
out_del_compressed_temp_plot=paste0(outfolder,in_template_name, "_del_data_for_compressed_temp_plot.csv")
out_ins_template_plot_data=paste0(outfolder,in_template_name, "_ins_template_plot_data.csv")
out_ins_template_plot_data_mechanism=paste0(outfolder,in_template_name, "_ins_template_plot_data_mechanism.csv")
out_consistency_true_ins=paste0(outfolder,in_template_name, "_SD-MMEJ_consistency_true_ins_del2.csv")


names(a)[names(a) == "ALIGNED_SEQ"] <- "RECONSTRUCTED_SEQ"
a10 <- subset(a, READS>9)
#Determine percent of total and inaccurate reads for each junction
a10$percent <- a10$READS/sum(a10$READS)*100
sum(a10$percent)
a10.nexact <- subset(a10, CLASS!="exact")
sum(a10.nexact$percent_inaccurate)
a10.exact <- subset(a10,CLASS=="exact")
a10.exact$percent_inaccurate <- "NA"
a10.nexact$percent_inaccurate <- a10.nexact$READS/sum(a10.nexact$READS)*100
a10 <- rbind(a10.exact,a10.nexact)
a10d<-subset(a10, CLASS_final=="deletion")
a10c<-subset(a10, CLASS_final=="complex")
a10i<-subset(a10, CLASS_final=="insertion")
a10e<-subset(a10, CLASS_final=="exact")
cc<-rbind(a10d,a10c,a10i,a10e)

write.csv(cc, out_cc) # VARIABLE

a10ci<-rbind(a10c, a10i)
del<-read.table(in_consistency_log, sep="",header = T) # File to be generated from the python output _consistency_log file input into excel as a comma/tab/space delimited file
# ERROR: more columns than column names

del$SEQ <- toupper(del$SEQ)
names(del)[names(del) == "SEQ"] <- "RECONSTRUCTED_SEQ"
dmerge<- merge(del, a10d,"RECONSTRUCTED_SEQ", all=TRUE)
library(tidyverse)
view(del)
write.csv(dmerge, out_dmerge) # VARIABLE

ins<-read.csv(in_ins_consistency) 
comp<-read.csv(in_comp_consistency) # VARIABLE


ic<-rbind(ins, comp)
icmerge<-merge(ic, a10ci, "RECONSTRUCTED_SEQ", all=FALSE)

# Combine deletion and insertion consistency tables to create single table that has all of the consistency information on it with accurate deletion boundaries.
ins1 <- icmerge[,19:34]
ins1$RECONSTRUCTED_SEQ <- icmerge$RECONSTRUCTED_SEQ
ins1$READS <- icmerge$READS
ins1$MICROHOMOLOGY <- icmerge$MICROHOMOLOGY 
ins1$MH_Length <- icmerge$MH_Length
ins1$NUMBER_OF_ALIGNMENTS  <- icmerge$NUMBER_OF_ALIGNMENTS
ins1$percent <- icmerge$percent
ins1$percent_inaccurate <- icmerge$percent_inaccurate
ins1$REPAIR_TYPE <- "InDel"
ins1$CLASS <- "Insertion"
ins1$left_del <- icmerge$left_del 
ins1$right_del <- icmerge$right_del
ins1$CONSISTENCY <- icmerge$consistency

del1 <- cbind(dmerge[,10:25])
del1$INSERTED_SEQ <- NA
del1$RECONSTRUCTED_SEQ <- dmerge$RECONSTRUCTED_SEQ
del1$READS <- dmerge$READS
del1$MICROHOMOLOGY <- dmerge$MICROHOMOLOGY 
del1$MH_Length <- dmerge$MH_Length
del1$NUMBER_OF_ALIGNMENTS  <- dmerge$NUMBER_OF_ALIGNMENTS
del1$percent <- dmerge$percent
del1$percent_inaccurate <- dmerge$percent_inaccurate
del1$REPAIR_TYPE <- dmerge$REPAIR_TYPE
del1 <- cbind(del1, dmerge[,6:8])
colnames(del1)
colnames(ins1)
names(del1) <- names(ins1)
all <- rbind(ins1,del1)
colnames(all)
# sum(all$percent_inaccurate)
write.csv(all,out_all) # VARIABLE


# FINAL FIGURES
# ====================================================== Insertion Primer Plot ====================================================================
# csv inputs - _combined_curated and _insertion_consistency
# csv outputs - _break

rm(list=ls())
library(stringr)
library(Biostrings)
library(ggplot2)
library(plyr)
library(grid)
#setwd("/Users/awelterofcollidingmaterials/desktop/Nick Woodward - Cas9_Iw7 Control Injections/6 Plots")

m <- out_cc #read.csv(out_cc) # VARIABLE
a <- in_ins_consistency #read.csv(in_ins_consistency) # VARIABLE
b <- in_comp_consistency #read.csv(in_comp_consistency) # VARIABLE

# colnames(a)
# colnames(b)
ab <-rbind(a,b)
# Input length of the full sequence
ab$RIGHT_DEL <- 127-(nchar(as.character(ab$RECONSTRUCTED_SEQ))-ab$right_del+1)
# 325-(nchar(as.character(ab$RECONSTRUCTED_SEQ)))
# ab$right_del+1

# Subset insertions based on plasmid and genotype
ins1 <- subset(m, CLASS_final=="insertion") # MAY WANT TO CHANGE TO CLASS_final
comp<-subset(m, CLASS_final=="complex")
ins<-rbind(ins1,comp)
ab.t <- subset(ab, consistency==TRUE)
colnames(ins)[19] <- "RECONSTRUCTED_SEQ"
mutIw7a <- merge(ins, ab.t, by="RECONSTRUCTED_SEQ")


# DIRECT REPEAT - RIGHT SIDE
mutIw7a.dr <- mutIw7a[!is.na(mutIw7a$DR_START),]
mutIw7a.dr$side <- ifelse(mutIw7a.dr$DR_END>mutIw7a.dr$right_del,
                          "RIGHT", "LEFT")
test <- as.data.frame(mutIw7a.dr[,42]) # the loop out column
test2 <- as.data.frame(sapply(test,gsub,pattern="-", replacement=""))
mutIw7a.dr<- cbind(as.data.frame(mutIw7a.dr),test2) 
colnames(mutIw7a.dr)[46] <- "RM"
mutIw7a.dr <- mutIw7a.dr[!is.na(mutIw7a.dr$RM),]
mutIw7a.dr$insertion_length <- nchar(as.character(mutIw7a.dr$insertion))
# Input length of the full sequence
mutIw7a.dr$DR_START_TRUE <- ifelse(mutIw7a.dr$side=="RIGHT", 
                                   (127-(nchar(as.character(mutIw7a.dr$RECONSTRUCTED_SEQ))-mutIw7a.dr$DR_START-mutIw7a.dr$insertion_length)),
                                   mutIw7a.dr$DR_START) # Gets the true start point for the direct repeat
# if right side it's 325-actual sequence length-DR_START-insertion length, if left its the DR_START
mutIw7a.dr$DR_END_TRUE <- ifelse(mutIw7a.dr$side=="RIGHT", 
                                 (127-(nchar(as.character(mutIw7a.dr$RECONSTRUCTED_SEQ))- mutIw7a.dr$DR_END-mutIw7a.dr$insertion_length)),
                                 mutIw7a.dr$DR_END) # Gets the true end point for the direct repeat
mutIw7a.dr$motif_ID <- paste(mutIw7a.dr$ID.y,mutIw7a.dr$DR_START_TRUE,mutIw7a.dr$DR_END_TRUE,sep="-")
mutIw7a.dr <- mutIw7a.dr[!duplicated(mutIw7a.dr$motif_ID),] # Removes all duplicate repeat motif IDs

drr <- subset(mutIw7a.dr, side=="RIGHT") # Subsetting the direct repeats from only the right side
a6 <- data.frame(ID=numeric() ,x=numeric() , y=numeric() , stringsAsFactors = FALSE) # Creates new data frame "a6"
for (i in 1:nrow(drr)){ 
  ins <- as.character(drr[i,38])    # creates insertion sequence region to search using "insertion" column
  dr1 <- as.data.frame(str_locate_all(drr[i,46], ins)) # search deletion sequence for string using "RM" column
  if (nrow(dr1)!=1){
    dr2 <- dr1
    while (nrow(dr2)!=1){
      for (j in 1:nrow(dr1)){
        if (dr1[j,2]==nchar(allowNA=FALSE, as.character(drr[i,46]))){
          break
        }else{
          prim <- paste(substring(drr[i,46],1,dr1[j,1]-1), substring(drr[i,46],dr1[j,2]+1,nchar(as.character(drr[i,46]))),sep="")
          s <- dr1[j,1]-1
          der <- substring(as.character(drr[i,37]), first=drr[i,35]-(s-1), last=drr[i,36]+(nchar(as.character(drr[i,46]))-dr1[j,2]))
          dr3 <- as.data.frame(str_locate_all(der, prim))
          if (nrow(dr3)==0){ # if there is a match, then make dr3=dr2, which will make the nrow(dr2)=1, breaking the while loop
            next
          } else{
            dr2 <- dr3
            k <- j
          }
        }
      }
      break
    }
    a7 <- cbind(drr[i,50], dr1[k,]) # once we know which one we want, we put that dr1 row
    a6 <- rbind(a6,a7)
  } else {
    a7 <- cbind(drr[i,50], dr1)
    a6 <- rbind(a6,a7)
  }
}

a6$is_trans <- ifelse(drr$DR_START_TRUE<77,
                      ifelse(drr$DR_START_TRUE>77, "trans", "no"), "no") 
a6$p1_start <- ifelse(a6$is_trans=="trans", # primer1 start need to be true (trans) if we're going to count it
                      drr$DR_START_TRUE, 
                      (drr$DR_END_TRUE+1)-(nchar(as.character(drr$RM))-a6[,3]))
a6$p1_end <- ifelse(a6$is_trans=="trans",
                    drr$DR_START_TRUE+(a6[,2]-2),
                    drr$DR_END_TRUE)
a6$p1_length <- 1+a6$p1_end-a6$p1_start
a6$mh_start <- ifelse(a6$is_trans=="trans",
                      (drr$DR_END_TRUE+1)-(nchar(as.character(drr$RM))-a6[,3]),
                      drr$DR_START_TRUE)
a6$mh_end <- ifelse(a6$is_trans=="trans",
                    drr$DR_END_TRUE,
                    drr$DR_START_TRUE+(a6[,2]-2))
a6$mh_length <- 1+a6$mh_end-a6$mh_start
a6$ins_start <- drr$DR_START_TRUE+(a6[,2]-1)
a6$ins_end <- (drr$DR_END_TRUE)-(nchar(as.character(drr$RM))-a6[,3])
a6$insertion_length <- 1+a6$ins_end-a6$ins_start 
a6$p2_start <- ifelse(a6$is_trans=="trans",
                      a6$p1_start,
                      drr$RIGHT_DEL+1)
a6$p2_end <-   ifelse(a6$is_trans=="trans",
                      a6$p1_end,
                      (drr$RIGHT_DEL+1)+(a6$p1_length-1))
a6$p2_length <- 1+a6$p2_end-a6$p2_start
colnames(a6)[1] <- "motif_ID"
drr2 <- merge(drr,a6, by="motif_ID")
# see if I can add the p2 that have the same p2 indexes
drr2$p2_ID <- paste(drr2$p2_start, drr2$p2_end, sep="-")
drr.u.ag <- aggregate(percent_inaccurate~p2_ID,data=drr2, sum)
drr.u.t <- as.data.frame(table(drr2$p2_ID))
drr.u.agt <- merge(drr.u.ag, drr.u.t, by.x="p2_ID", by.y="Var1")
drr.um <- merge(drr.u.agt,drr2, by="p2_ID")
drr.umu <- drr.um[!duplicated(drr.um$p2_ID),]
drr.umu$p2_start2 <- drr.umu$p2_start-.5
drr.umu$p2_end2 <- drr.umu$p2_end+.5
p <- ggplot(drr.umu, aes(x=p2_ID))
p + geom_boxplot(aes(ymin = p2_start2, lower = p2_start2, middle = p2_start2, 
                     upper = p2_end2, ymax = p2_end2,
                     colour=percent_inaccurate.x, fill=percent_inaccurate.x), 
                 stat = "identity") + 
  coord_flip()+
  scale_colour_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                         values=c(1.0,0.8,0.6,0.4,0.2,0)) +
  scale_fill_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                       values=c(1.0,0.8,0.6,0.4,0.2,0))+
  scale_y_continuous(limits=c(130,190),
                     breaks=51:125,
                     labels = c("A", "G", "C", "C", "G","A", "A", "T", "T", "C", "G", "G", "T", "A", "C", "A", "T", "T", "A", "C", "C", 
                                "C", "T", "G", "t", "t", "a", "t", "C", "C", "C", "T", "A", "G", "C", 
                                "G", "G", "C", "C", "G", "C", "A","T", "A","G", "G", "C", "C", 
                                "A", "C", "T","A", "G", "T", "G", "G", "A", "T","C", "T", "G", "G", "A", 
                                "T", "C", "C","T","C","T","A","G","A","G","T","C"), #CCCTAGC  TATGGTC   TGCGCTACT   AGTGGATCTGGGGCC   GCATAGGCC
                     name = "Nucleotide",
                     expand=c(0,0))+
  theme_classic(base_size = 18)+
  theme(axis.text.y= element_blank())+
  geom_hline(yintercept = 74.5, colour="red")+
  geom_hline(yintercept = 78.5,colour="red")

# DIRECT REPEAT - LEFT SIDE
drl <- subset(mutIw7a.dr, mutIw7a.dr$side=="LEFT")
a6 <- data.frame(ID=numeric() ,x=numeric() , y=numeric() , stringsAsFactors = FALSE)
for (i in 1:nrow(drl)){
  ins <- as.character(drl[i,38])    # create insertion sequence region to search
  dr1 <- as.data.frame(str_locate_all(drl[i,46], ins)) # search deletion sequence for string
  if (nrow(dr1)!=1){
    dr2 <- dr1
    while (nrow(dr2)!=1){
      for (j in 1:nrow(dr1)){
        if (dr1[j,2]==nchar(allowNA = TRUE,as.character(drl[i,46]))){
          break
        }else{
          prim <- paste(substring(drl[i,46],1,dr1[j,1]-1), substring(drl[i,46],dr1[j,2]+1,nchar(as.character(drl[i,46]))),sep="")
          s <- dr1[j,1]-1
          der <- substring(as.character(drl[i,37]), first=drl[i,35]-(s-1), last=drl[i,36]+(nchar(as.character(drl[i,46]))-dr1[j,2]))
          dr3 <- as.data.frame(str_locate_all(der, prim))
          if (nrow(dr3)==0){ # if there is a match, then make dr3=dr2, which will make the nrow(dr2)=1, breaking the while loop
            next
          } else{
            dr2 <- dr3
            k <- j
          }
        }
      }
      break
    }
    a7 <- cbind(drl[i,50], dr1[k,]) # once we know which one we want, we put that dr1 row
    a6 <- rbind(a6,a7)
  } else {
    a7 <- cbind(drl[i,50], dr1)
    a6 <- rbind(a6,a7)
  }
}
a6$is_trans <- ifelse(drl$DR_END_TRUE<161,
                      ifelse(drl$DR_END_TRUE>161, "trans", "no"), "no")
a6$p1_start <- ifelse(a6$is_trans=="trans",
                      (drl$DR_END_TRUE+1)-(nchar(as.character(drl$RM))-a6[,3]),
                      drl$DR_START_TRUE)
a6$p1_end <- ifelse(a6$is_trans=="trans",
                    drl$DR_END_TRUE,
                    drl$DR_START_TRUE+(a6[,2]-2))
a6$p1_length <- 1+a6$p1_end-a6$p1_start
a6$mh_start <- ifelse(a6$is_trans=="trans",
                      drl$DR_START_TRUE,
                      (drl$DR_END_TRUE+1)-(nchar(as.character(drl$RM))-a6[,3]))
a6$mh_end <- ifelse(a6$is_trans=="trans",
                    drl$DR_START_TRUE+(a6[,2]-2),
                    drl$DR_END_TRUE)
a6$mh_length <- 1+a6$mh_end-a6$mh_start
a6$ins_start <- drl$DR_START_TRUE+(a6[,2]-1)
a6$ins_end <- (drl$DR_END_TRUE)-(nchar(as.character(drl$RM))-a6[,3])
a6$insertion_length <- 1+a6$ins_end-a6$ins_start 
a6$p2_start <- ifelse(a6$is_trans=="trans",
                      a6$p1_start,
                      (drl$left_del)-(a6$p1_length-1))
a6$p2_end <-   ifelse(a6$is_trans=="trans",
                      a6$p1_end,
                      drl$left_del)
a6$p2_length <- 1+a6$p2_end-a6$p2_start
colnames(a6)[1] <- "motif_ID"
drl2 <- merge(drl,a6, by="motif_ID")
drl2$p2_ID <- paste(drl2$p2_start, drl2$p2_end, sep="-")
drl2.u.ag <- aggregate(percent_inaccurate~p2_ID,data=drl2, sum)
drl2.u.t <- as.data.frame(table(drl2$p2_ID))
drl2.u.agt <- merge(drl2.u.ag, drl2.u.t, by.x="p2_ID", by.y="Var1")
drl2.um <- merge( drl2.u.agt,drl2, by="p2_ID")
drl2.umu <- drl2.um[!duplicated(drl2.um$p2_ID),]
drl2.umu$p2_start2 <- drl2.umu$p2_start-.5
drl2.umu$p2_end2 <- drl2.umu$p2_end+.5
p <- ggplot(drl2.umu, aes(x=p2_ID))
p + geom_boxplot(aes(ymin = p2_start2, lower = p2_start2, middle = p2_start2, 
                     upper = p2_end2, ymax = p2_end2,
                     colour=percent_inaccurate.x, fill=percent_inaccurate.x), 
                 stat = "identity") + 
  coord_flip()+
  scale_colour_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                         values=c(1.0,0.8,0.6,0.4,0.2,0)) +
  scale_fill_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                       values=c(1.0,0.8,0.6,0.4,0.2,0))+
  scale_y_continuous(limits=c(130,190),
                     breaks=51:125,
                     labels = c("A", "G", "C", "C", "G","A", "A", "T", "T", "C", "G", "G", "T", "A", "C", "A", "T", "T", "A", "C", "C", 
                                "C", "T", "G", "t", "t", "a", "t", "C", "C", "C", "T", "A", "G", "C", 
                                "G", "G", "C", "C", "G", "C", "A","T", "A","G", "G", "C", "C", 
                                "A", "C", "T","A", "G", "T", "G", "G", "A", "T","C", "T", "G", "G", "A", 
                                "T", "C", "C","T","C","T","A","G","A","G","T","C"), #CCCTAGC  TATGGTC   TGCGCTACT   AGTGGATCTGGGGCC   GCATAGGCC
                     name = "Nucleotide",
                     expand=c(0,0))+
  theme_classic(base_size = 18)+
  theme(axis.text.y= element_blank())+
  geom_hline(yintercept = 75, colour="red")+
  geom_hline(yintercept = 78,colour="red")


# REVERSE COMPLEMENT - LEFT SIDE
mutIw7a.rc <- mutIw7a[!is.na(mutIw7a$RC_START),]
mutIw7a.rc$side <- ifelse(mutIw7a.rc$RC_END>mutIw7a.rc$right_del,
                          "RIGHT", "LEFT")
test <- as.data.frame(mutIw7a.rc[,43]) # snap back column
test2 <- as.data.frame(sapply(test,gsub,pattern="-", replacement=""))
mutIw7a.rc <- cbind(as.data.frame(mutIw7a.rc),test2)
colnames(mutIw7a.rc)[46] <- "RM"
mutIw7a.rc$insertion_length <- nchar(as.character(mutIw7a.rc$insertion))
mutIw7a.rc$RC_START_TRUE <- ifelse(mutIw7a.rc$side=="RIGHT", 
                                   (127-(nchar(as.character(mutIw7a.rc$RECONSTRUCTED_SEQ))- mutIw7a.rc$RC_START-mutIw7a.rc$insertion_length)),
                                   mutIw7a.rc$RC_START)
mutIw7a.rc$RC_END_TRUE <- ifelse(mutIw7a.rc$side=="RIGHT", 
                                 (127-(nchar(as.character(mutIw7a.rc$RECONSTRUCTED_SEQ))- mutIw7a.rc$RC_END-mutIw7a.rc$insertion_length)),
                                 mutIw7a.rc$RC_END)
mutIw7a.rc$motif_ID <- paste(mutIw7a.rc$ID.y,mutIw7a.rc$RC_START_TRUE,mutIw7a.rc$RC_END_TRUE,sep="-")
mutIw7a.rc <- mutIw7a.rc[!duplicated(mutIw7a.rc$motif_ID),]
rcl <- subset(mutIw7a.rc, side=="LEFT")

a6 <- data.frame(ID=numeric() ,x=numeric() , y=numeric() , stringsAsFactors = FALSE)
for (i in 1:nrow(rcl)){
  ins <- as.character(rcl[i,38])    # create insertion sequence region to search
  dnains <- lapply(ins, DNAString)
  revins <- lapply(dnains, reverseComplement)
  ins <- lapply(revins, toString)
  rm <- toupper(rcl[i,46])
  dr1 <- as.data.frame(str_locate_all(rm, as.character(ins))) # search deletion sequence for string
  if (nrow(dr1)!=1){
    dr2 <- dr1
    while (nrow(dr2)!=1){
      for (j in 1:nrow(dr1)){
        if (dr1[j,2]==nchar(as.character(rcl[i,46]))){
          break
        }else{
          prim <- paste(substring(rcl[i,46],1,dr1[j,1]-1), substring(rcl[i,46],dr1[j,2]+1,nchar(as.character(rcl[i,46]))),sep="")
          dnains2 <- lapply(prim, DNAString)
          revins2 <- lapply(dnains2, reverseComplement)
          prim <- lapply(revins2, toString)
          s <- dr1[j,1]-1
          der <- substring(as.character(rcl[i,37]), first=rcl[i,35]-(s-1), last=rcl[i,36]+(nchar(as.character(rcl[i,46]))-dr1[j,2]))
          dr3 <- as.data.frame(str_locate_all(der, as.character(prim)))
          if (nrow(dr3)==0){ # if there is a match, then make dr3=dr2, which will make the nrow(dr2)=1, breaking the while loop
            next
          } else{
            dr2 <- dr3
            k <- j
          }
        }
      }
      break
    }
    a7 <- cbind(rcl[i,50], dr1[k,]) # once we know which one we want, we put that dr1 row
    a6 <- rbind(a6,a7)
  } else {
    a7 <- cbind(rcl[i,50], dr1)
    a6 <- rbind(a6,a7)
  }
}
a6$is_trans <- "no"
a6$p1_start <- (rcl$RC_END_TRUE)-(nchar(as.character(rcl$RM))-a6[,3]-1)
a6$p1_end <- rcl$RC_END_TRUE
a6$p1_length <- 1+a6$p1_end-a6$p1_start
a6$mh_start <- rcl$RC_START_TRUE
a6$mh_end <- rcl$RC_START_TRUE+(a6[,2]-2)
a6$mh_length <- 1+a6$mh_end-a6$mh_start
a6$ins_start <- rcl$RC_START_TRUE+(a6[,2]-1)
a6$ins_end <- (rcl$RC_START_TRUE)+(a6[,3]-1)
a6$insertion_length <- 1+a6$ins_end-a6$ins_start 
a6$p2_start <- (rcl$left_del)-(a6$p1_length-1)
a6$p2_end <- rcl$left_del
a6$p2_length <- 1+a6$p2_end-a6$p2_start
colnames(a6)[1] <- "motif_ID"
rcl2 <- merge(rcl,a6, by="motif_ID")
rcl2$p2_ID <- paste(rcl2$p2_start, rcl2$p2_end, sep="-")
rcl.u.ag <- aggregate(READS~p2_ID,data=rcl2, sum)
rcl.u.t <- as.data.frame(table(rcl2$p2_ID))
rcl.u.agt <- merge(rcl.u.ag, rcl.u.t, by.x="p2_ID", by.y="Var1")
rcl.um <- merge( rcl.u.agt,rcl2, by="p2_ID")
rcl.umu <- rcl.um[!duplicated(rcl.um$p2_ID),]
rcl.umu$p2_start2 <- rcl.umu$p2_start-.5
rcl.umu$p2_end2 <- rcl.umu$p2_end+.5


# REVERSE COMPLEMENT - RIGHT SIDE
rcr <- subset(mutIw7a.rc, side=="RIGHT")
a6 <- data.frame(ID=numeric() ,x=numeric() , y=numeric() , stringsAsFactors = FALSE)
for (i in 1:nrow(rcr)){
  ins <- as.character(rcr[i,38])    # create insertion sequence region to search
  dnains <- lapply(ins, DNAString)
  revins <- lapply(dnains, reverseComplement)
  ins <- lapply(revins, toString)
  rm <- toupper(rcr[i,46])
  dr1 <- as.data.frame(str_locate_all(rm, as.character(ins))) # search deletion sequence for string
  if (nrow(dr1)!=1){
    dr2 <- dr1
    while (nrow(dr2)!=1){
      for (j in 1:nrow(dr1)){
        if (dr1[j,2]==nchar(as.character(rcr[i,46]))){
          break
        }else{
          prim <- paste(substring(rcr[i,46],1,dr1[j,1]-1), substring(rcr[i,46],dr1[j,2]+1,nchar(as.character(rcr[i,46]))),sep="")
          dnains2 <- lapply(prim, DNAString)
          revins2 <- lapply(dnains2, reverseComplement)
          prim <- lapply(revins2, toString)
          s <- dr1[j,1]-1
          der <- substring(as.character(rcr[i,37]), first=rcr[i,35]-(s-1), last=rcr[i,36]+(nchar(as.character(rcr[i,46]))-dr1[j,2]))
          dr3 <- as.data.frame(str_locate_all(der, as.character(prim)))
          if (nrow(dr3)==0){ # if there is a match, then make dr3=dr2, which will make the nrow(dr2)=1, breaking the while loop
            next
          } else{
            dr2 <- dr3
            k <- j
          }
        }
      }
      break
    }
    a7 <- cbind(rcr[i,50], dr1[k,]) # once we know which one we want, we put that dr1 row
    a6 <- rbind(a6,a7)
  } else {
    a7 <- cbind(rcr[i,50], dr1)
    a6 <- rbind(a6,a7)
  }
}
a6$is_trans <- "no"
a6$p1_start <- rcr$RC_START_TRUE
a6$p1_end <- rcr$RC_START_TRUE+(a6[,2]-2)
a6$p1_length <- 1+a6$p1_end-a6$p1_start
a6$mh_start <-  (rcr$RC_END_TRUE)-(nchar(as.character(rcr$RM))-a6[,3]-1)
a6$mh_end <- rcr$RC_END_TRUE 
a6$mh_length <- 1+a6$mh_end-a6$mh_start
a6$ins_start <- rcr$RC_START_TRUE+(a6[,2]-1)
a6$ins_end <- (rcr$RC_START_TRUE)+(a6[,3]-1)
a6$insertion_length <- 1+a6$ins_end-a6$ins_start 
a6$p2_start <- rcr$RIGHT_DEL+1
a6$p2_end <- (rcr$RIGHT_DEL+1)+(a6$p1_length-1)
a6$p2_length <- 1+a6$p2_end-a6$p2_start
colnames(a6)[1] <- "motif_ID"
rcr2 <- merge(rcr,a6, by="motif_ID")
# see if I can add the p2 that have the same p2 indexes
rcr2$p2_ID <- paste(rcr2$p2_start, rcr2$p2_end, sep="-")
rcr.u.ag <- aggregate(READS~p2_ID,data=rcr2, sum)
rcr.u.t <- as.data.frame(table(rcr2$p2_ID))
rcr.u.agt <- merge(rcr.u.ag, rcr.u.t, by.x="p2_ID", by.y="Var1")
rcr.um <- merge( rcr.u.agt,rcr2, by="p2_ID")
rcr.umu <- rcr.um[!duplicated(rcr.um$p2_ID),]
rcr.umu$p2_start2 <- rcr.umu$p2_start-.5
rcr.umu$p2_end2 <- rcr.umu$p2_end+.5
p <- ggplot(rcr.umu, aes(x=p2_ID))
p + geom_boxplot(aes(ymin = p2_start2, lower = p2_start2, middle = p2_start2, 
                     upper = p2_end2, ymax = p2_end2,
                     colour=percent_inaccurate, fill=percent_inaccurate), 
                 stat = "identity") + 
  coord_flip()+
  scale_colour_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                         values=c(1.0,0.8,0.6,0.4,0.2,0)) +
  scale_fill_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                       values=c(1.0,0.8,0.6,0.4,0.2,0))+
  scale_y_continuous(limits=c(130,190),
                     breaks=51:125,
                     labels = c("A", "G", "C", "C", "G","A", "A", "T", "T", "C", "G", "G", "T", "A", "C", "A", "T", "T", "A", "C", "C", 
                                "C", "T", "G", "t", "t", "a", "t", "C", "C", "C", "T", "A", "G", "C", 
                                "G", "G", "C", "C", "G", "C", "A","T", "A","G", "G", "C", "C", 
                                "A", "C", "T","A", "G", "T", "G", "G", "A", "T","C", "T", "G", "G", "A", 
                                "T", "C", "C","T","C","T","A","G","A","G","T","C"), #CCCTAGC  TATGGTC   TGCGCTACT   AGTGGATCTGGGGCC   GCATAGGCC
                     name = "Nucleotide",
                     expand=c(0,0))+
  theme_classic(base_size = 18)+
  theme(axis.text.y= element_blank())+
  geom_hline(yintercept = 75, colour="red")+
  geom_hline(yintercept = 78,colour="red")

mutIw7dr <- rbind(drr2,drl2)
mutIw7dr1 <- mutIw7dr[,c(1:2,5:50)]
colnames(mutIw7dr1)
mutIw7dr1$RC_START_TRUE <- "NA"
mutIw7dr1$RC_END_TRUE <- "NA"
mutIw7dr2 <- cbind(mutIw7dr1,mutIw7dr[,51:66])
mutIw7dr2$mechanism <- ifelse(mutIw7dr2$is_trans=="trans",
                              "Trans",
                              "Loop-out")
mutIw7rc <- rbind(rcr2,rcl2)
mutIw7rc1 <- mutIw7rc[,c(1:2,5:48)]
mutIw7rc1$DR_START_TRUE <- "NA"
mutIw7rc1$DR_END_TRUE <- "NA"
mutIw7rc1<-cbind(mutIw7rc1,mutIw7rc$RC_START_TRUE)
mutIw7rc1<-cbind(mutIw7rc1,mutIw7rc$RC_END_TRUE)
mutIw7rc2 <- cbind(mutIw7rc1,mutIw7rc[,51:66])
mutIw7rc2$mechanism <- "Snap-back"
names(mutIw7rc2) <- names(mutIw7dr2)
mutIw7break <- rbind(mutIw7rc2, mutIw7dr2)
# mutIw7break<-mutIw7dr2
colnames(mutIw7rc2)
colnames(mutIw7dr2)
mutIw7break$p2_mechID <- paste(mutIw7break$p2_ID,mutIw7break$mechanism,sep="-")
mutIw7.u.ag <- aggregate(READS~p2_mechID,data=mutIw7break, sum)
mutIw7.u.t <- as.data.frame(table(mutIw7break$p2_mechID))
mutIw7.u.agt <- merge(mutIw7.u.ag, mutIw7.u.t, by.x="p2_mechID", by.y="Var1")
mutIw7.um <- merge( mutIw7.u.agt,mutIw7break, by="p2_mechID")
mutIw7.umu <- mutIw7.um[!duplicated(mutIw7.um$p2_mechID),]

mutIw7.umu$p2_start2 <- mutIw7.umu$p2_start-.5
mutIw7.umu$p2_end2 <- mutIw7.umu$p2_end+.5
mutIw7.umu$p1_start2 <- mutIw7.umu$p1_start-.5
mutIw7.umu$p1_end2 <- mutIw7.umu$p1_end+.5
mutIw7.umu$mh_start2 <- mutIw7.umu$mh_start-.5
mutIw7.umu$mh_end2 <- mutIw7.umu$mh_end+.5
mutIw7.umu$ins_start2 <- mutIw7.umu$ins_start-.5
mutIw7.umu$ins_end2 <- mutIw7.umu$ins_end+.5
mutIw7.umu$percent_insertion_jxn <- mutIw7.umu$READS.x/sum(mutIw7.umu$READS.x)*100
f2 <- ggplot()+
  geom_boxplot(data=mutIw7.umu, aes(x=p2_mechID, ymin = p2_start2, lower = p2_start2, middle = p2_start2, 
                                    upper = p2_end2, ymax = p2_end2,
                                    fill=percent_insertion_jxn, linetype=mechanism),
               colour="black",
               stat = "identity",
               width=.8) + 
  coord_flip()+
  scale_linetype_manual(values=c("solid", "dashed","dotted"),guide=FALSE)+
  scale_fill_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                       values=c(1.0,0.8,0.6,0.4,0.2,0),
                       limits=c(0,16),
                       guide_legend(title="% Single-step\nInsertion Reads"))+
  scale_y_continuous(limits=c(130,190),
                     breaks=51:125,
                     labels = c("A", "G", "C", "C", "G","A", "A", "T", "T", "C", "G", "G", "T", "A", "C", "A", "T", "T", "A", "C", "C", 
                                "C", "T", "G", "t", "t", "a", "t", "C", "C", "C", "T", "A", "G", "C", 
                                "G", "G", "C", "C", "G", "C", "A","T", "A","G", "G", "C", "C", 
                                "A", "C", "T","A", "G", "T", "G", "G", "A", "T","C", "T", "G", "G", "A", 
                                "T", "C", "C","T","C","T","A","G","A","G","T","C"), #CCCTAGC  TATGGTC   TGCGCTACT   AGTGGATCTGGGGCC   GCATAGGCC
                     name = "Nucleotide",
                     expand=c(0,0))+
  theme_classic(base_size = 18)+
  theme(axis.text.y= element_blank(),
        axis.title.y=element_blank(),
        axis.ticks=element_blank())+
  geom_rect(data=rectangles, aes(ymin=as.numeric(xmin), ymax=as.numeric(xmax), xmin=ymin, xmax=ymax),
            fill="grey", alpha=0.4, colour=NA)+
  geom_boxplot(data=mutIw7.umu, aes(x=p2_mechID, ymin = p2_start2, lower = p2_start2, middle = p2_start2, 
                                    upper = p2_end2, ymax = p2_end2,
                                    fill=percent_insertion_jxn, linetype=mechanism),
               colour="white",
               stat = "identity",
               width=.8)+  
  geom_hline(yintercept = c(74.5,78.5), colour="red", size=0.75)
f2

f2 <- ggplot()+
  geom_boxplot(data=mutIw7.umu, aes(x=p2_mechID, ymin = p2_start2, lower = p2_start2, middle = p2_start2, 
                                    upper = p2_end2, ymax = p1_start,
                                    fill=percent_insertion_jxn, linetype=mechanism),
               colour="black",
               stat = "identity",
               width=.8) + 
  coord_flip()+
  scale_linetype_manual(values=c("solid", "dashed","dotted"),guide=FALSE)+
  scale_fill_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                       values=c(1.0,0.8,0.6,0.4,0.2,0),
                       guide_legend(title="% Inaccurate\nReads"))+
  scale_y_continuous(limits=c(130,190),
                     breaks=51:125,
                     labels = c("A", "G", "C", "C", "G","A", "A", "T", "T", "C", "G", "G", "T", "A", "C", "A", "T", "T", "A", "C", "C", 
                                "C", "T", "G", "t", "t", "a", "t", "C", "C", "C", "T", "A", "G", "C", 
                                "T", "A", "T", "G", "G", "T", "C", 
                                "T", "G", "C", "G", "C", "T", "A", "C", "T", 
                                "A", "G", "T", "G", "G", "A", "T", "C", "T", "G", "G", "G", "G", "C", "C",
                                "G","C","A","T","A","G","G","C","C"), #CCCTAGC  TATGGTC   TGCGCTACT   AGTGGATCTGGGGCC   GCATAGGCC
                     name = "Nucleotide",
                     expand=c(0,0))+
  theme_classic(base_size = 18)+
  theme(axis.text.y= element_blank(),
        axis.title.y=element_blank(),
        axis.ticks=element_blank())+
  geom_rect(data=rectangles, aes(ymin=as.numeric(xmin), ymax=as.numeric(xmax), xmin=ymin, xmax=ymax),
            fill="grey", alpha=0.4, colour=NA)+
  geom_boxplot(data=mutIw7.umu, aes(x=p2_mechID, ymin = p2_start2, lower = p2_start2, middle = p2_start2, 
                                    upper = p2_end2, ymax = p2_end2,
                                    fill=percent_insertion_jxn, linetype=mechanism),
               colour="white",
               stat = "identity",
               width=.8)+ 
  geom_boxplot(data=mutIw7.umu, aes(x=p2_mechID, ymin = p1_start2, lower = p1_start2, middle = p1_start2, 
                                    upper = p1_end2, ymax = p1_end2,
                                    fill=percent_insertion_jxn, linetype=mechanism),
               colour="white",
               stat = "identity",
               width=.8)+
  geom_boxplot(data=mutIw7.umu, aes(x=p2_mechID, ymin = ins_start2, lower = ins_start2, middle = ins_start2, 
                                    upper = ins_end2, ymax = ins_end2,
                                    fill=percent_insertion_jxn, linetype=mechanism),
               colour="white",
               stat = "identity",
               width=.8)+
  geom_hline(yintercept = c(74.5,78.5), colour="red", size=0.8)
f2


# Create _break.csv
write.csv (mutIw7break, out_break)
mutIw7break<-read.csv(in_break_curated)
attach(mutIw7break)

#create primer set unique indices to aggregate over
mutIw7break$jxnID <- paste(p1_start,p1_end,p2_start,p2_end, mechanism,sep="-")
mutIw7jxn.ag <- aggregate(READS~jxnID,data=mutIw7break, sum)
mutIw7jxn.t <- as.data.frame(table(mutIw7break$jxnID))
mutIw7jxn.agt <- merge(mutIw7jxn.ag, mutIw7jxn.t, by.x="jxnID", by.y="Var1")
mutIw7jxn.um <- merge( mutIw7jxn.agt,mutIw7break, by="jxnID")
mutIw7jxn.umu <- mutIw7jxn.um[!duplicated(mutIw7jxn.um$jxnID),]
mutIw7jxn.umu$p2_start2 <- mutIw7jxn.umu$p2_start-.5
mutIw7jxn.umu$p2_end2 <- mutIw7jxn.umu$p2_end+.5
mutIw7jxn.umu$p1_start2 <- mutIw7jxn.umu$p1_start-.5
mutIw7jxn.umu$p1_end2 <- mutIw7jxn.umu$p1_end+.5
mutIw7jxn.umu$mh_start2 <- mutIw7jxn.umu$mh_start-.5
mutIw7jxn.umu$mh_end2 <- mutIw7jxn.umu$mh_end+.5
mutIw7jxn.umu$ins_start2 <- mutIw7jxn.umu$ins_start-.5
mutIw7jxn.umu$ins_end2 <- mutIw7jxn.umu$ins_end+.5
mutIw7jxn.umu$percent_insertion_jxn <- mutIw7jxn.umu$READS.x/sum(mutIw7jxn.umu$READS.x)*100
mutIw7jxn.uml <- subset(mutIw7jxn.umu, side=="LEFT")
mutIw7jxn.umr <- subset(mutIw7jxn.umu, side=="RIGHT")
mutIw7jxn.uml <- mutIw7jxn.uml[order(mutIw7jxn.uml$p2_start),]
mutIw7jxn.umr <- mutIw7jxn.umr[order(mutIw7jxn.umr$p2_start, decreasing=TRUE),]
mutIw7jxn.uml$order <- paste(10:(9+nrow(mutIw7jxn.uml)), "-", sep="")   
mutIw7jxn.umr$order <- paste(10:(9+nrow(mutIw7jxn.umr)), "-", sep="")   # if the number of rows is greater than the number of letters, then dont put [1:nrow(...)]
mutIw7jxn.um2 <- rbind(mutIw7jxn.uml, mutIw7jxn.umr)
mutIw7jxn.um2 <- mutIw7jxn.um2[order(mutIw7jxn.um2$order,decreasing=TRUE),]

mutIw7jxn.um2 <- arrange(transform(mutIw7jxn.um2,
                                   jxnID=factor(jxnID,levels=jxnID)),jxnID)

#make a pretty plot
rect_left <- c(33.5,43.5,53.5,63.5,73.5,83.5,93.5,103.5,113.5,123.5,133.5,143.5)
rectangles <- data.frame(
  xmin = rect_left,
  xmax = rect_left + 5,
  ymin = -Inf,
  ymax = Inf
)
rectangles
rectangles[9,2] <- 118.5

f2 <- ggplot()+
  geom_boxplot(data=mutIw7jxn.um2, aes(x=jxnID, ymin = p2_start2, lower = p2_start2, middle = p2_start2, 
                                       upper = p2_end2, ymax = p1_start,
                                       fill=percent_insertion_jxn, linetype=mechanism),
               colour="black",
               stat = "identity",
               width=.8) + 
  coord_flip()+
  scale_linetype_manual(values=c("solid", "dashed","dotted"),guide=FALSE)+
  scale_fill_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                       values=c(1.0,0.8,0.6,0.4,0.2,0),
                       limits=c(0,13),
                       guide_legend(title="% Single-step\nInsertion Reads"))+
  scale_y_continuous(limits=c(40,120),
                     breaks=44:115,
                     labels = c("A", "G", "C", "A", #static
                                "G", "C", "C", "G", "A", "A", "T", "T", "C", "G", "G", "T", "A",
                                "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                                "C", "G", "G", "C", "C", "G", "C", "A", "T", "A", "G", "G", "C", "C",
                                "A", "C", "T", "A", "G", "T", "G", "G", "A", "T", "C", "T", "G",
                                "G", "A", "T", "C", "C", "T", "C", "T"),
                    
                     # Iw7 AGCA
                     # AGCA
                     # GCCGAATTCGGTA
                     # CATTACCCTGTTATCCCTAG
                     # CGGCCGCATAGGCC
                     # ACTAGTGGATCTG
                     # GATCCTCT
                     
                     # GCAGCCGAATT
                     # CGGTACATTACCCTGTTATCCCTAG
                     # CGGCCGCATAGGCC
                     # ACTAGTGGATCTGGAT
                     # FLSBMH
                     # labels = c("A", "G", "C", "A", #static
                     #            "T", "C", "A", "C", "T", "A", "C", "G", "C", "T", "T", "A", "C",
                     #            "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                     #            "G", "G", "G", "C", "C", "T", "A", "T", "A", "T", "G", "G", "C", "C",
                     #            "G", "C", "C", "C", "G", "T", "C", "A", "T", "C", "C", "C", "T",
                     #            "G", "A", "T", "C", "C", "T", "C", "T"),
                     # 
                     # PolyT1 
                     # labels = c("A", "G", "C", "A", #static
                     #            "C", "T", "A", "T", "A", "C", "A", "A", "C", "A", "T", "C", "G",
                     #            "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                     #            "C", "G", "G", "C", "C", "T", "T", "T", "T", "T", "G", "G", "C", "C",
                     #            "G", "G", "T", "T", "C", "A", "T", "C", "T", "G", "G", "A", "A",
                     #            "G", "A", "T", "C", "C", "T", "C", "T"),
                     # PolyT2
                     # labels = c("A", "G", "C", "A", #static
                     #            "T", "A", "G", "A", "C", "C", "C", "T", "A", "T", "T", "G", "T",
                     #            "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                     #            "A", "G", "G", "C", "C", "T", "T", "T", "T", "T", "G", "G", "C", "C",
                     #            "C", "C", "T", "A", "T", "C", "A", "A", "T", "A", "A", "T", "T",
                     #            "G", "A", "T", "C", "C", "T", "C", "T"),
                     # PolyT3 
                     # labels = c("A", "G", "C", "A", #static
                     #            "C", "G", "A", "C", "A", "A", "G", "G", "C", "G", "C", "A", "G",
                     #            "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                     #            "A", "G", "G", "C", "C", "T", "T", "T", "T", "T", "G", "G", "C", "C",
                     #            "G", "G", "T", "T", "T", "C", "T", "A", "A", "A", "T", "A", "A",
                     #            "G", "A", "T", "C", "C", "T", "C", "T"),
                     # PolyT4 
                     # labels = c("A", "G", "C", "A", #static
                     #            "G", "C", "G", "A", "G", "T", "A", "A", "T", "T", "T", "T", "G",
                     #            "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                     #            "A", "G", "G", "C", "C", "T", "T", "T", "T", "T", "G", "G", "C", "C",
                     #            "C", "G", "C", "T", "C", "A", "A", "C", "G", "G", "T", "A", "G",
                     #            "G", "A", "T", "C", "C", "T", "C", "T"),
                     
                     # LO2SB3
                     # labels = c("A", "G", "C", "A", #static
                     #            "G", "A", "G", "A", "G", "A", "T", "T", "C", "G", "C", "A", "T",
                     #            "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                     #            "A", "G", "G", "C", "C", "T", "C", "G", "C", "C", "G", "G", "C", "C",
                     #            "G", "G", "G", "C", "T", "C", "C", "T", "T", "C", "A", "T", "T",
                     #            "G", "A", "T", "C", "C", "T", "C", "T"),
                     
                     # LO5P2SB2
                     # labels = c("A", "G", "C", "A", #static
                     #            "G", "G", "G", "G", "T", "C", "T", "T", "A", "G", "C", "C", "T",
                     #            "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                     #            "T", "G", "G", "C", "C", "G", "T", "A", "C", "C", "G", "G", "C", "C",
                     #            "G", "G", "A", "G", "A", "G", "C", "C", "C", "G", "A", "T", "T",
                     #            "G", "A", "T", "C", "C", "T", "C", "T"),
                     
                     # LO2   
                     # labels = c("A", "G", "C", "A", #static
                     #            "A", "A", "G", "A", "T", "T", "A", "G", "G", "C", "C", "T", "G",
                     #            "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                     #            "C", "G", "G", "C", "C", "T", "C", "A", "C", "C", "G", "G", "C", "C",
                     #            "C", "G", "C", "G", "G", "A", "G", "C", "G", "A", "A", "T", "C",
                     #            "G", "A", "T", "C", "C", "T", "C", "T"),
                     
                     # LO3
                     # labels = c("A", "G", "C", "A", #static
                     #            "A", "T", "A", "T", "G", "G", "T", "C", "G", "G", "T", "T", "C",
                     #            "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                     #            "T", "G", "G", "C", "C", "T", "A", "C", "C", "C", "G", "G", "C", "C",
                     #            "C", "T", "C", "C", "A", "T", "A", "T", "A", "A", "T", "T", "T",
                     #            "G", "A", "T", "C", "C", "T", "C", "T"),
                     
                     # LO5
                     # labels = c("A", "G", "C", "A", #static
                     #            "T", "C", "T", "T", "T", "T", "A", "C", "A", "A", "C", "T", "A",
                     #            "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                     #            "A", "G", "G", "C", "C", "A", "T", "C", "C", "C", "G", "G", "C", "C",
                     #            "T", "T", "G", "C", "C", "G", "C", "C", "T", "T", "A", "T", "A",
                     #            "G", "A", "T", "C", "C", "T", "C", "T"),
                     
                     # # LO2_4DEL
                     # labels = c("A", "G", "C", "A", #static
                     #            "C", "G", "T", "A", "G", "T", "T", "A", "A", "G", "A", "T", "C",
                     #            "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                     #            "T", "G", "G", "C", "C", "A", "A", "T", "C", "C", "G", "G", "C", "C",
                     #            "T", "A", "A", "T", "T", "T", "C", "G", "C", "T", "T", "C", "C",
                     #            "G", "A", "T", "C", "C", "T", "C", "T"),
                     # 
                     # SB2
                     # labels = c("A", "G", "C", "A", #static
                     #            "C", "T", "C", "A", "C", "C", "T", "T", "A", "T", "T", "G", "T",
                     #            "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                     #            "A", "G", "G", "C", "C", "A", "A", "A", "C", "A", "G", "G", "C", "C",
                     #            "G", "G", "C", "G", "C", "C", "T", "C", "C", "T", "A", "A", "C",
                     #            "G", "A", "T", "C", "C", "T", "C", "T"),
                     
                     # SB4
                     # labels = c("A", "G", "C", "A", #static
                     #            "G", "A", "A", "T", "C", "A", "G", "T", "C", "G", "A", "T", "T",
                     #            "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                     #            "A", "G", "G", "C", "C", "A", "A", "T", "A", "T", "G", "G", "C", "C",
                     #            "G", "G", "G", "A", "A", "T", "C", "T", "A", "A", "C", "C", "A",
                     #            "G", "A", "T", "C", "C", "T", "C", "T"),
                     
# AGCATTACCCTGTTATCCC-TAGAGGCCTTTTTGGCCGGTTTCTAAATAAGATCCTCT
                     name = "Nucleotide",
                     expand=c(0,0))+
  theme_classic(base_size = 10)+
  theme(axis.text.y= element_blank(),
        axis.title.y=element_blank(),
        axis.ticks=element_blank())+
  geom_rect(data=rectangles, aes(ymin=as.numeric(xmin), ymax=as.numeric(xmax), xmin=ymin, xmax=ymax),
            fill="grey", alpha=0.4, colour=NA)+
  geom_boxplot(data=mutIw7jxn.umu, aes(x=jxnID, ymin = p2_start2, lower = p2_start2, middle = p2_start2, 
                                       upper = p2_end2, ymax = p2_end2,
                                       fill=percent_insertion_jxn, linetype=mechanism),
               colour="black",
               stat = "identity",
               width=.8)+ 
  geom_boxplot(data=mutIw7jxn.um2, aes(x=jxnID, ymin = p1_start2, lower = p1_start2, middle = p1_start2, 
                                       upper = p1_end2, ymax = p1_end2,
                                       fill=percent_insertion_jxn, linetype=mechanism),
               colour="black",
               stat = "identity",
               width=.8)+
  geom_hline(yintercept = 77.5, colour="red", size=0.6)
f2


# ================================================ Insertion Repeat Motif Plot (Expanded) ==================================================================
# csv input - _break
# csv output - _insertion_consistent_analysis
rm(list=ls())
library(plyr)
library(grid)
library(ggplot2)
library(dplyr)
library(Biostrings)
library(stringr)
library(dplyr)
#setwd("/Users/awelterofcollidingmaterials/desktop/Nick Woodward - Cas9_Iw7 Control Injections/6 Plots")

mutIw7break <- out_break #read.csv(out_break)
mutIw7break$MOTIF_START <- ifelse(mutIw7break$DR_START_TRUE=="NA",
                                  mutIw7break$RC_START_TRUE,
                                  mutIw7break$DR_START_TRUE)
mutIw7break$MOTIF_END <- ifelse(mutIw7break$DR_END_TRUE=="NA",
                                mutIw7break$RC_END_TRUE,
                                mutIw7break$DR_END_TRUE)
mutIw7break$motif <- paste(mutIw7break$MOTIF_START,mutIw7break$MOTIF_END,sep="-")
mutIw7break$motif_mechID <- paste(mutIw7break$motif,mutIw7break$mechanism,sep="-")

#save it
write.csv(mutIw7break,out_ins_consistency)
# mutIw7break <- read.csv("LO5P2SB2_insertion_consistent_analysis.csv")

mutIw7break<-read.csv(out_ins_consistency)
im.ag <- aggregate(READS~motif_ID,data=mutIw7break,sum)
im.t <- as.data.frame(table(mutIw7break$motif_ID))
im.agt <- merge(im.ag,im.t,by.x="motif_ID", by.y="Var1")
im.um <- merge(mutIw7break,im.agt,by="motif_ID")
im.uml <- subset(im.um, side=="LEFT")
im.umr <- subset(im.um,side=="RIGHT")
im.uml <- im.uml[order(as.numeric(im.uml$MOTIF_START),decreasing=TRUE),]
im.umr <- im.umr[order(as.numeric(im.umr$MOTIF_START),decreasing=TRUE),]
im.uml$order <- paste(10:(9+nrow(im.uml)), "-", sep="")
im.umr$order <- paste(10:(9+nrow(im.umr)), "-", sep="")   # if the number of rows is greater than the number of letters, then dont put [1:nrow(...)]
im.um2 <- rbind(im.uml, im.umr)
im.um2 <- im.um2[order(im.um2$order,decreasing=TRUE),]
im.um3 <- im.um2[!duplicated(im.um2$motif_mechID),]

neworder <- im.um3$motif_ID
im.um4 <- arrange(transform(im.um3,
                            motif_ID=factor(motif_ID,levels=motif_ID)),motif_ID)
# im.um4<-im.um3
im.um4$percent_insertion <- im.um4$READS.y/sum(im.um4$READS.y)*100
im.um4$MOTIF_START2 <- as.numeric(im.um4$MOTIF_START)-.5
im.um4$MOTIF_END2 <- as.numeric(im.um4$MOTIF_END)+.5

rect_left <- c(33.5,43.5,53.5,63.5,73.5,83.5,93.5,103.5,113.5,123.5,133.5,143.5)
rectangles <- data.frame(
  xmin = rect_left,
  xmax = rect_left + 5,
  ymin = -Inf,
  ymax = Inf
)
rectangles
rectangles[9,2] <- 118.5

f2 <- ggplot()+
  geom_boxplot(data=im.um4, aes(x=motif_ID, ymin = MOTIF_START2, lower = MOTIF_START2, middle = MOTIF_START2, 
                                upper = MOTIF_END2, ymax = MOTIF_END2,
                                fill=percent_insertion, linetype=mechanism),
               colour="white",
               stat = "identity",
               width=.8,
               lwd=.5) + 
  coord_flip()+
  scale_linetype_manual(values=c("solid", "dashed","dotted"),guide=FALSE)+
  scale_fill_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                       values=c(1.0,0.8,0.6,0.4,0.2,0),
                       guide_legend(title="% Single-step\nInsertion Reads"),
                       limits=c(0,9.25))+
  scale_y_continuous(limits=c(40,120),
                     breaks=44:115,
                     labels = c("A", "G", "C", "A", #static
                                "G", "C", "C", "G", "A", "A", "T", "T", "C", "G", "G", "T", "A",
                                "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                                "C", "G", "G", "C", "C", "G", "C", "A", "T", "A", "G", "G", "C", "C",
                                "A", "C", "T", "A", "G", "T", "G", "G", "A", "T", "C", "T", "G",
                                "G", "A", "T", "C", "C", "T", "C", "T"),
                     name = "Nucleotide",
                     expand=c(0,0))+
  theme_classic(base_size = 10)+
  theme(legend.text=element_text(size=8,face="bold"),
        axis.text.y= element_blank(),
        axis.title.y=element_blank(),
        axis.text.x= element_text(face="bold"),
        axis.title.x=element_text(face="bold"),
        axis.ticks=element_blank(),
        legend.position=c(0,1), legend.justification=c(.05,2.31),
        legend.key.size=unit(0.22,"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  guides(guide_legend(title.theme = element_text(size=8, face="bold", angle=0),
                      label.theme=element_text(size=8, face="bold", angle=0)))+
  geom_rect(data=rectangles, aes(ymin=as.numeric(xmin), ymax=as.numeric(xmax), xmin=ymin, xmax=ymax),
            fill="grey", alpha=0.4, colour=NA)+
  geom_boxplot(data=im.um4, aes(x=motif_ID, ymin = MOTIF_START2, lower = MOTIF_START2, middle = MOTIF_START2, 
                                upper = MOTIF_END2, ymax = MOTIF_END2,
                                fill=percent_insertion, linetype=mechanism),
               colour="black",
               stat = "identity",
               width=.8,
               lwd=.5)+  
  geom_hline(yintercept = 77.5, colour="red", size=0.75)
f2
# dev.off()
# save as "mut_m4_ins_rm_plot_08102017.pdf" 7.00787x4.33071 inches





# ============================================= Deletion Repeat Motif Plot (Expanded) ======================================================
# csv input - _combined_curated, _python_spreadsheet (data from top of python output), _sample_ID (ID and RECONSTRUCTED_SEQ of python input), del_data_for_template_plot
# csv output - del_data_for_template_plot, del_data_for_compressed_temp_plot-mech, del_data_for_compressed_temp_plot, 
rm(list=ls())
library(plyr)
library(grid)
library(ggplot2)
# setwd("/Users/awelterofcollidingmaterials/Desktop/Iw7R_HiSeq/Flexible Loop Analysis")
m <- read.csv(out_cc)
# colnames(m)[18] <- "RECONSTRUCTED_SEQ"
del <- subset(m, CLASS=="deletion")
a <- read.csv("Iw7_python_spreadsheet.csv") # This file is from the deletion program's python consistency table. Import into excel from txt and then take the top of it.
b <- read.csv("Iw7_sample_ID.csv") # this file is created by opening the python output file and copying the sequences and ID into an excel sheet and splitting it into ID and Reconstructued seq --> to get reconstructed seq, remove all "-" in the aligned seqs

ab <- merge(a, b, by.x="Sample.ID", by.y="ID")
abc <- merge(ab,m, by="RECONSTRUCTED_SEQ", all.x=TRUE)

abc$start <- ifelse(abc$Break.Side=="left",
                    (78-abc$Motif.to.Break), (78+abc$Motif.to.Break-abc$Motif.Length))
abc$end <- ifelse(abc$Break.Side=="left",
                  (77-abc$Motif.to.Break+abc$Motif.Length), (77+abc$Motif.to.Break))
abc$motif_mechID <- paste(abc$start, abc$end, abc$Mechanism, sep="-")

abc.ag <- aggregate(READS~motif_mechID,data=abc,sum)
abc.t <- as.data.frame(table(abc$motif_mechID))
abc.agt <- merge(abc.ag,abc.t,by.x="motif_mechID", by.y="Var1")
abc.um <- merge(abc,abc.agt,by="motif_mechID")
abc.um <- abc.um[!duplicated(abc.um$motif_mechID),]
abc.uml <- subset(abc.um,start<77)
abc.umr <- subset(abc.um,start>77)
abc.uml <- abc.uml[order(abc.uml$start),]
abc.umr <- abc.umr[order(abc.umr$start),]
abc.uml$order <- paste(10:(9+nrow(abc.uml)), "-", sep="")   
abc.umr$order <- paste(10:(9+nrow(abc.umr)), "-", sep="")   # if the number of rows is greater than the number of letters, then dont put [1:nrow(...)]
abc.um2 <- rbind(abc.uml, abc.umr)
abc.um2 <- abc.um2[order(abc.um2$order,decreasing=TRUE),]

neworder <- abc.um2$motif_mechID
abc.um3 <- arrange(transform(abc.um2,
                             motif_mechID=factor(motif_mechID,levels=motif_mechID)),motif_mechID)
abc.um3$percent_deletion <- abc.um3$READS.y/sum(abc.um3$READS.y)*100
abc.um3$MOTIF_START2 <- as.numeric(abc.um3$start)-.5
abc.um3$MOTIF_END2 <- as.numeric(abc.um3$end)+.5

write.csv(abc.um3,out_del_for_template_plot)

sweet <- abc.um3
a6 <- data.frame(x=numeric() , y=numeric(), z=numeric(),w=numeric() , stringsAsFactors = FALSE)
for(i in 1:nrow(sweet)){
  g <- as.data.frame(sweet[i,50]:sweet[i,51])
  a7 <- cbind(sweet[i,1], g)
  a7 <- cbind(a7, sweet[i,46])
  a7 <- cbind(a7, sweet[i,6])
  a6 <- rbind(a6,a7)
}

colnames(a6) <- c("motif_mechID", "temp_coord","READS","mechanism" )
abcct.sb <- subset(a6, mechanism=="snap-back")
abcct.lo <- subset(a6, mechanism=="loop-out")
abcct.sb.ag <- aggregate(READS~temp_coord,data=abcct.sb, sum)
abcct.lo.ag <- aggregate(READS~temp_coord,data=abcct.lo, sum)
abcct.sb.ag$mechanism <- "Snap-back"
abcct.lo.ag$mechanism <- "Loop-out"
abcct.ag <- rbind(abcct.sb.ag, abcct.lo.ag)
abcct.ag$x <- "1"
abcct.ag$percent_deletion <- abcct.ag$READS/sum(abcct.ag$READS)*100
abcct.ag$temp_coord2 <- abcct.ag$temp_coord-0.5
abcct.ag2 <- aggregate(READS~temp_coord,data=a6, sum)
abcct.ag2$x <- "Iw7"
abcct.ag2$percent_deletion <- abcct.ag2$READS/sum(abcct.ag2$READS)*100
abcct.ag2$temp_coord2 <- abcct.ag2$temp_coord-0.5
summary(abcct.ag)
summary(abcct.ag2)
sum(abcct.ag$percent_deletion)

write.csv(abcct.ag,out_del_compressed_temp_plot_mech)
write.csv(abcct.ag2,out_del_compressed_temp_plot)

# rm(list=ls())
mutm4temp <- read.csv(out_del_for_template_plot)

# For the expanded template plot, i need rearrange stuff and also not include samples that are outside of my sequence range
mutm4temp <- arrange(transform(mutm4temp,
                           motif_mechID=factor(motif_mechID,levels=motif_mechID)),motif_mechID)
mutm4temp2 <- subset(mutm4temp, MOTIF_START2<143&MOTIF_END2>33)

rect_left <- c(33.5,43.5,53.5,63.5,73.5,83.5,93.5,103.5,113.5,123.5,133.5,143.5)
rectangles <- data.frame(
  xmin = rect_left,
  xmax = rect_left + 5,
  ymin = -Inf,
  ymax = Inf
)
rectangles
rectangles[9,2] <- 118.5
# Need to remove those jucntions that start less than 53....sadly
mutm4temp3 <- subset(mutm4temp2, MOTIF_END2>33)
mutm4temp4 <- subset(mutm4temp3, MOTIF_START2<143)

f2 <- ggplot()+
  geom_boxplot(data=mutm4temp4, aes(x=motif_mechID, ymin = MOTIF_START2, lower = MOTIF_START2, middle = MOTIF_START2, 
                                upper = MOTIF_END2, ymax = MOTIF_END2,
                                fill=percent_deletion, linetype=Mechanism),
               colour="white",
               stat = "identity",
               width=.8,  
               lwd=.5) + 
  coord_flip()+
  scale_linetype_manual(values=c("solid", "dashed","dotted"),guide=FALSE)+
  scale_fill_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                       values=c(1.0,0.8,0.6,0.4,0.2,0),
                       limits=c(0,14.5),
                       guide_legend(title="% SD-MMEJ\nConsistent\nDeletion Reads"))+
  scale_y_continuous(limits=c(40,120),
                     breaks=44:115,
                     labels = c("A", "G", "C", "A", #static
                                "G", "C", "C", "G", "A", "A", "T", "T", "C", "G", "G", "T", "A",
                                "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                                "C", "G", "G", "C", "C", "G", "C", "A", "T", "A", "G", "G", "C", "C",
                                "A", "C", "T", "A", "G", "T", "G", "G", "A", "T", "C", "T", "G",
                                "G", "A", "T", "C", "C", "T", "C", "T"),
                     name = "Nucleotide",
                     expand=c(0,0))+
  theme_classic(base_size = 10)+
  theme(legend.text=element_text(size=8,face="bold"),
        axis.text.y= element_blank(),
        axis.title.y=element_blank(),
        axis.text.x= element_text(face="bold"),
        axis.title.x=element_text(face="bold"),
        axis.ticks=element_blank(),
        legend.position=c(.1,.56), legend.justification=c(.4,2),
        legend.key.size=unit(0.42,"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  guides(guide_legend(title.theme = element_text(size=4, face="bold", angle=0),
                     label.theme=element_text(size=8, face="bold", angle=0)))+
  geom_rect(data=rectangles, aes(ymin=as.numeric(xmin), ymax=as.numeric(xmax), xmin=ymin, xmax=ymax),
            fill="grey", alpha=0.8, colour=NA)+
  geom_boxplot(data=mutm4temp4, aes(x=motif_mechID, ymin = MOTIF_START2, lower = MOTIF_START2, middle = MOTIF_START2, 
                                upper = MOTIF_END2, ymax = MOTIF_END2,
                                fill=percent_deletion, linetype=Mechanism),
               colour="black",
               stat = "identity",
               width=.8, 
               lwd=.5)+  
  geom_hline(yintercept = 77.5, colour="red", size=0.75)
f2
 #saved as "mut_m4_del_rm_plot-08082017.pdf" 4.33x7.00787






# ============================================== Insertion Repeat Motif Plots (Compressed) ==============================================================
rm(list=ls())
#setwd("/Users/awelterofcollidingmaterials/Desktop/sdmmej-master/Test_Data/PolyT1_output")
# csv input - _insertion_consistent_analysis, _ins_template_plot_data.csv
# csv output - _ins_template_plot_data.csv, ins_template_plot_data_mechanism.csv

# Need to sum the percent of inaccurate reads for EVERY nucleotide. So based off of the repeat motif indeces, make a line where every nucleotide of that motif has the same percent of inaccurate reads
im.um3<-read.csv(out_ins_consistency)
a6 <- data.frame(x=numeric() , y=numeric(), z=numeric(),w=numeric() , stringsAsFactors = FALSE)
for(i in 1:nrow(im.um3)){
  g <- as.data.frame(im.um3[i,71]:im.um3[i,72])
  a7 <- cbind(im.um3[i,3], g)
  a7 <- cbind(a7, im.um3[i,21])
  a7 <- cbind(a7, im.um3[i,69])
  a6 <- rbind(a6,a7)
}
colnames(a6) <- c("motif_ID", "temp_coord","READS","mechanism" )
mutm4ct.sb <- subset(a6, mechanism=="Snap-back")
mutm4ct.lo <- subset(a6, mechanism=="Loop-out")
# mutm4ct.tr <- subset(a6,mechanism=="Trans")
mutm4ct.sb.ag <- aggregate(READS~temp_coord,data=mutm4ct.sb, sum)
mutm4ct.lo.ag <- aggregate(READS~temp_coord,data=mutm4ct.lo, sum)
# mutm4ct.tr.ag <- aggregate(X._OF_READS~temp_coord,data=mutm4ct.tr, sum)
mutm4ct.sb.ag$mechanism <- "Snap-back"
mutm4ct.lo.ag$mechanism <- "Loop-out"
# mutm4ct.tr.ag$mechanism <- "Trans"
mutm4ct.ag <- rbind(mutm4ct.sb.ag, mutm4ct.lo.ag)
mutm4ct.ag$x <- "1"
mutm4ct.ag$percent_insertion <- mutm4ct.ag$READS/sum(mutm4ct.ag$READS)*100
mutm4ct.ag2 <- aggregate(READS~temp_coord,data=a6, sum)
mutm4ct.ag2$x <- "Iw7"
mutm4ct.ag2$percent_insertion <- mutm4ct.ag2$READS/sum(mutm4ct.ag2$READS)*100
summary(mutm4ct.ag)

write.csv(mutm4ct.ag, out_ins_template_plot_data, row.names=FALSE)
# added coordinates that were missing with 0 reads so that it will work correctly in the plot

mutm4ct.ag3<-mutm4ct.ag # read.csv(out_ins_template_plot_data)
mutm4ct.ag3$x <- in_template_name
mutm4ct.ag3$temp_coord2 <- mutm4ct.ag3$temp_coord-0.5

write.csv(mutm4ct.ag3,out_ins_template_plot_data_mechanism)

#make grey rectangles to put behind the plots
# rect_left <- c(133.5,143.5,153.5,163.5,173.5,183.5,193.5)
# rectangles <- data.frame(
#   xmin = rect_left,
#   xmax = rect_left + 5,
#   ymin = -Inf,
#   ymax = Inf
# )
# rectangles
# rectangles[7,2] <- 198.5

f1 <- ggplot() + 
  geom_line(data=mutm4ct.ag3,aes(x=as.numeric(temp_coord2), y=mechanism,colour=percent_insertion), size=25)+
  scale_colour_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                         values=c(1.0,0.8,0.6,0.4,0.2,0),
                         guide_legend(title="% Single-step\nInsertion Reads"),
                         limits=c(0,6.5))+
  scale_x_continuous(limits=c(40,120),
                     breaks=44:115,
                     labels = c("A", "G", "C", "A", #static
                                "G", "C", "C", "G", "A", "A", "T", "T", "C", "G", "G", "T", "A",
                                "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                                "C", "G", "G", "C", "C", "G", "C", "A", "T", "A", "G", "G", "C", "C",
                                "A", "C", "T", "A", "G", "T", "G", "G", "A", "T", "C", "T", "G",
                                "G", "A", "T", "C", "C", "T", "C", "T"),
                     # PolyT1 
                     # labels = c("A", "G", "C", "A", #static
                     #            "C", "T", "A", "T", "A", "C", "A", "A", "C", "A", "T", "C", "G",
                     #            "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                     #            "C", "G", "G", "C", "C", "T", "T", "T", "T", "T", "G", "G", "C", "C",
                     #            "G", "G", "T", "T", "C", "A", "T", "C", "T", "G", "G", "A", "A",
                     #            "G", "A", "T", "C", "C", "T", "C", "T"),
                     # PolyT2
                     # labels = c("A", "G", "C", "A", #static
                     #            "T", "A", "G", "A", "C", "C", "C", "T", "A", "T", "T", "G", "T",
                     #            "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                     #            "A", "G", "G", "C", "C", "T", "T", "T", "T", "T", "G", "G", "C", "C",
                     #            "C", "C", "T", "A", "T", "C", "A", "A", "T", "A", "A", "T", "T",
                     #            "G", "A", "T", "C", "C", "T", "C", "T"),
                     # PolyT3 
                     # labels = c("A", "G", "C", "A", #static
                     #            "C", "G", "A", "C", "A", "A", "G", "G", "C", "G", "C", "A", "G",
                     #            "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                     #            "A", "G", "G", "C", "C", "T", "T", "T", "T", "T", "G", "G", "C", "C",
                     #            "G", "G", "T", "T", "T", "C", "T", "A", "A", "A", "T", "A", "A",
                     #            "G", "A", "T", "C", "C", "T", "C", "T"),
                     # PolyT4 
                     # labels = c("A", "G", "C", "A", #static
                     #            "G", "C", "G", "A", "G", "T", "A", "A", "T", "T", "T", "T", "G",
                     #            "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                     #            "A", "G", "G", "C", "C", "T", "T", "T", "T", "T", "G", "G", "C", "C",
                     #            "C", "G", "C", "T", "C", "A", "A", "C", "G", "G", "T", "A", "G",
                     #            "G", "A", "T", "C", "C", "T", "C", "T"),
                     
                     expand=c(0,0))+
  theme_classic(base_size = 9)+
  theme(legend.text=element_text(face="bold"),
        strip.text.x=element_text(face="bold"),
        axis.title.y=element_blank(),
        axis.text.y=element_text(face="bold"),
        axis.text.x= element_text(face="bold"),
        axis.title.x=element_blank(),
        axis.ticks=element_blank(),
        legend.key.size=unit(0.25,"cm"),
        plot.margin = unit(c(.1,0,0,0), "cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  geom_line(data=mutm4ct.ag3,aes(x=as.numeric(temp_coord2), y=mechanism,colour=percent_insertion), size=4)+
  geom_vline(xintercept = 77.1, colour="red", size=0.75)
f1



#=================================================== Deletion Repeat Motif Plot (Compressed) ============================================================
# IMPORTANT: make sure you go through and insert the missing coordinates and give them 0 reads and percent_deletions. Or else the graph will be wrong
# csv input - _del_data_for_compressed_temp_plot-mech, 
# csv output - none
rm(list=ls())
#setwd("/Users/awelterofcollidingmaterials/Desktop/sdmmej-master/Test_Data/PolyT4_output")

mutm4ctm <- read.csv(out_del_compressed_temp_plot_mech)
mutm4ctm$temp_coord2 <- mutm4ctm$temp_coord2+0.5
f1 <- ggplot() + 
  geom_line(data=mutm4ctm,aes(x=as.numeric(temp_coord2), y=mechanism,colour=percent_deletion), size=20)+
  scale_colour_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                         values=c(1.0,0.8,0.6,0.4,0.2,0),
                         guide_legend(title="% SD-MMEJ\nConsistent\nDeletion Reads"),
                         limits=c(0,3.5))+
  scale_x_continuous(limits=c(40,120),
                     breaks=44:115,
                     labels = c("A", "G", "C", "A", #static
                                "G", "C", "C", "G", "A", "A", "T", "T", "C", "G", "G", "T", "A",
                                "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                                "C", "G", "G", "C", "C", "G", "C", "A", "T", "A", "G", "G", "C", "C",
                                "A", "C", "T", "A", "G", "T", "G", "G", "A", "T", "C", "T", "G",
                                "G", "A", "T", "C", "C", "T", "C", "T"),
                     # PolyT1 
                     # labels = c("A", "G", "C", "A", #static
                     #            "C", "T", "A", "T", "A", "C", "A", "A", "C", "A", "T", "C", "G",
                     #            "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                     #            "C", "G", "G", "C", "C", "T", "T", "T", "T", "T", "G", "G", "C", "C",
                     #            "G", "G", "T", "T", "C", "A", "T", "C", "T", "G", "G", "A", "A",
                     #            "G", "A", "T", "C", "C", "T", "C", "T"),
                     # PolyT2
                     # labels = c("A", "G", "C", "A", #static
                     #            "T", "A", "G", "A", "C", "C", "C", "T", "A", "T", "T", "G", "T",
                     #            "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                     #            "A", "G", "G", "C", "C", "T", "T", "T", "T", "T", "G", "G", "C", "C",
                     #            "C", "C", "T", "A", "T", "C", "A", "A", "T", "A", "A", "T", "T",
                     #            "G", "A", "T", "C", "C", "T", "C", "T"),
                     # PolyT3 
                     # labels = c("A", "G", "C", "A", #static
                     #            "C", "G", "A", "C", "A", "A", "G", "G", "C", "G", "C", "A", "G",
                     #            "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                     #            "A", "G", "G", "C", "C", "T", "T", "T", "T", "T", "G", "G", "C", "C",
                     #            "G", "G", "T", "T", "T", "C", "T", "A", "A", "A", "T", "A", "A",
                     #            "G", "A", "T", "C", "C", "T", "C", "T"),
                     # PolyT4 
                     # labels = c("A", "G", "C", "A", #static
                     #            "G", "C", "G", "A", "G", "T", "A", "A", "T", "T", "T", "T", "G",
                     #            "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                     #            "A", "G", "G", "C", "C", "T", "T", "T", "T", "T", "G", "G", "C", "C",
                     #            "C", "G", "C", "T", "C", "A", "A", "C", "G", "G", "T", "A", "G",
                     #            "G", "A", "T", "C", "C", "T", "C", "T"),
                     expand=c(0,0))+
  theme_classic(base_size = 10)+
  theme(legend.text=element_text(face="bold"),
        strip.text.x=element_text(face="bold"),
        axis.title.y=element_blank(),
        axis.text.y=element_text(face="bold"),
        axis.text.x= element_text(face="bold"),
        axis.title.x=element_blank(),
        axis.ticks=element_blank(),
        legend.key.size=unit(0.25,"cm"),
        plot.margin = unit(c(.1,0,0,0), "cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  geom_vline(xintercept = c(161.5), colour="red", size=0.75)
f1


# ========================================================= Deletion Stopper Plots ===================================================================
# csv input - _all_SD-MMEJ_consistency (analysis_annotated), _python_consistency (original python consistency file), _del_merge_table (From analysis annotated)
# csv output - _SD-MMEJ_consistency_true_ins_del2 (with total deletion and insertion true)
rm(list=ls())
#setwd("/Users/awelterofcollidingmaterials/desktop/Nick Woodward - Cas9_Iw7 Control Injections/6 Plots")
library(plyr)
library(grid)
library(ggplot2)
a <- read.csv(out_all)
a$PG <- paste("Iw7") 
ins <- subset(a,REPAIR_TYPE=="InDel")
a$TOTAL_DELETION_TRUE <- a$right_del-a$left_del
a$TOTAL_DELETION_TRUE <- ifelse(a$TOTAL_DELETION_TRUE>=0,
                                a$TOTAL_DELETION_TRUE,0)
a$TOTAL_INSERTION_TRUE <- a$READ_LENGTH+a$TOTAL_DELETION_TRUE-127
a$TOTAL_INSERTION_TRUE <- ifelse(a$REPAIR_TYPE=="MHJ",
                                 0,a$TOTAL_INSERTION_TRUE)
#save this
write.csv(a, out_consistency_true_ins)

a <- read.csv(in_consistency_log)
names(a)[names(a) == 'SEQ'] <- 'RECONSTRUCTED_SEQ'
a$RECONSTRUCTED_SEQ <- toupper(a$RECONSTRUCTED_SEQ)

PolyA1.merge <- read.csv(out_dmerge)

testa <- read.csv(out_all)
testa$RECONSTRUCTED_SEQ <- toupper(testa$RECONSTRUCTED_SEQ)
testa$PG <- "Iw7"

PolyA1 <- subset(a, !duplicated(a$RECONSTRUCTED_SEQ))
names(PolyA1.merge)[names(PolyA1.merge) == 'ALIGNED_SEQ'] <- 'RECONSTRUCTED_SEQ'
PolyA1.merge2 <- merge(PolyA1.merge, a, by="RECONSTRUCTED_SEQ", all=TRUE)
PolyA1merged <- merge(PolyA1.merge2, testa, by="RECONSTRUCTED_SEQ", all=FALSE)

# Take InDel data and make left and right boundaries, things that are further than 96 or 63bp will be a + on the graph
ins <- subset(testa, CLASS=="Insertion")

# Right side
ins.stop <- subset(ins, right_del<95)
ins.xtra <- subset(ins,right_del>95)
ins.xtra$right_del <- 95
ins.stop <- rbind(ins.stop,ins.xtra) # So far only formatted all of the deletion boundaries that will show up as + on the plots
# Now for the left side
ins.l.stop <- subset(ins.stop, left_del>57) 
ins.l.extras <- subset(ins.stop, left_del<57) # Less than 63 resolves as + on the left side of the deletion boundary graph
summary(ins.l.extras)
ins.l.extras$left_del <- 57
summary(ins.l.extras)
ins.stop <- rbind(ins.l.stop, ins.l.extras)
summary(ins.stop)

del<-read.csv(out_consistency_true_ins)

del<-cbind(PolyA1merged)
del.stop <- subset(del, RIGHT_DEL_INDEX.x<95)
del.xtra <- subset(del, RIGHT_DEL_INDEX.x>95)
summary(del.xtra)
del.xtra$RIGHT_DEL_INDEX.x <- 95
del.stop <- rbind(del.stop,del.xtra) # Resolved all right deletion indices outside of the graph limits (106) to be + on the graph
del.l.stop <- subset(del.stop, LEFT_DEL_INDEX.x>57)
del.l.extras <- subset(del.stop, LEFT_DEL_INDEX.x<57)
summary(del.l.extras)
del.l.extras$LEFT_DEL_INDEX.x <- 57
summary(del.l.extras)
del.stop <- rbind(del.l.stop, del.l.extras)
summary(del.stop)
pretty.del <- as.data.frame(del.stop[,64])
pretty.del$right_del <- del.stop[,65]
pretty.del$REPAIR <- "deletion"
pretty.del$CONSISTENCY <- del.stop[,66]
pretty.del$percent <- del.stop[,62]
pretty.del$PG <- del.stop$PG
colnames(pretty.del)[1:2] <- c("left_del","right_del")

#make Indel table same as pretty.del so we can combine and aggrigate
pretty.ins <- as.data.frame(ins.stop$left_del)
pretty.ins$right_del <- ins.stop$right_del
pretty.ins$REPAIR <- "InDel"
pretty.ins$CONSISTENCY <- ins.stop$CONSISTENCY
pretty.ins$percent <- ins.stop$percent_inaccurate
pretty.ins$PG <- ins.stop$PG
colnames(pretty.ins)[1:2] <- c("left_del","right_del")

all.del <- rbind(pretty.ins, pretty.del)
all.del.ls<-sort(all.del$left_del, decreasing = FALSE)
all.del.rs<-sort(all.del$right_del, decreasing = TRUE)
#Combine InDel and deletion data AND aggregate
del.agl <- aggregate(percent~left_del+CONSISTENCY+PG+REPAIR,data=all.del, sum)
del.agl$side <- "Left"
del.agr <- aggregate(percent~right_del+CONSISTENCY+PG+REPAIR,data=all.del, sum)
del.agr$side <- "Right"
colnames(del.agr)[1] <- "del_bound"
colnames(del.agl)[1] <- "del_bound"
del.ag <- rbind(del.agl, del.agr)

del.ag.false <- subset(del.ag, CONSISTENCY=="FALSE")
del.ag.true <- subset(del.ag, CONSISTENCY=="TRUE")
del.order <- rbind(del.ag.true, del.ag.false)

# rm(list=ls())
# setwd("/Users/awelterofcollidingmaterials/Desktop/sdmmej-master/Test_Data/PolyT4_output")
write.csv(del.order,"Iw7delorder.csv") #write csv file and manually change L and R flank depending on deletion boundary. Add percent inaccurate to get a nice plot
del.order<-read.csv("Iw7delorder.csv")

# Right side
right <- subset(del.order, side=="Right")
PolyA1db <- subset(right, PG=="Iw7")

f2 <- ggplot(PolyA1db, aes(x=del_bound, y=percent, alpha = CONSISTENCY, fill=REPAIR))+
  scale_y_continuous(name = "% Deletion Boundaries",
                     expand=c(0,0),
                     breaks=c(0,10,20,30,40,50),
                     limits=c(0,35))+
  scale_fill_brewer(palette="Dark2")+
  scale_alpha_discrete(range=c(0.2,1))+
  guides(fill=guide_legend(title="Repair\nOutcome"))+
  theme_classic(base_size = 18)+ 
  theme(axis.text.x=element_text(size=18, face="bold"),
        axis.text.y=element_text(size=18, face="bold"),
        axis.title.x=element_text(size=18, face="bold"))+
  theme(strip.text.x = element_text(size=20, face="bold"),
        strip.text.y = element_text(size=20, face="bold"))+
  geom_bar(stat="identity", colour="black") +
  scale_x_continuous(limits=c(40,120),
                     breaks=44:115,
                     labels = c("A", "G", "C", "A", #static
                                "G", "C", "C", "G", "A", "A", "T", "T", "C", "G", "G", "T", "A",
                                "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C", "T", "A", "G", #static
                                "C", "G", "G", "C", "C", "G", "C", "A", "T", "A", "G", "G", "C", "C",
                                "A", "C", "T", "A", "G", "T", "G", "G", "A", "T", "C", "T", "G",
                                "G", "A", "T", "C", "C", "T", "C", "T"),
                     # PolyT1 
                     # scale_x_continuous(limits=c(161,184),
                     #                    breaks=162:184,
                     #                    labels = c("T","A","G","C","G","G","C","C","T","T","T","T"
                     #                               ,"T","G","G","C","C","G","G","T","T","C","A"),
                     # PolyT3 
                     # labels = c("T", "A", "G","A", "G", "G", "C", "C", "T", "T", "T", "T", "T", "G", "G", "C", "C", "G", "G", "T", "T", "T", "C"),
                     
                     # PolyT4 
                     # labels = c("T", "A", "G", A", "G", "G", "C", "C", "T", "T", "T", "T", "T", "G", "G", "C", "C", "C", "G", "C", "T", "C", "A"),
                     
                     # AT1
                     # scale_x_continuous(limits=c(161,178),
                     #                    breaks=162:184,
                     #                    labels = c("T", "A", "G", "T", "G", "G", "C", "C", "A", "T", "A", "T", "A", "G", "G", "C", "C", "C", "T", "A", "C", "G", "A"),
                     # AT2
                     #scale_x_continuous(limits=c(161,178),
                     # breaks=162:184,
                     # labels = c("T", "A", "G", "G", "G", "G", "C", "C", "A", "T", "A", "T", "A", "G", "G", "C", "C", "A", "T", "T", "G", "T", "A"),
                     # AT3 TAGTGGCCATATAGGCCCGAACT
                     # AT4 TAGTGGCCTATATGGCCTGTTCA
                     # AT5 TAGGGGCCTATATGGCCGCCCGT
                     # AT6 TAGTGGCCTATATGGCCTCTATC
                     # AT7 TAGGGGCCTATATGGCCCTCTATGGCC
                     
                     expand=c(0,0))
f2

# Left side
left <- subset(del.order, side=="Left")
mutIw7CRdbl <- subset(left, PG=="AT7")

f2 <- ggplot(mutIw7CRdbl, aes(x=del_bound, y=percent, alpha = CONSISTENCY, fill=REPAIR))+
  scale_y_continuous(name = "% Deletion Boundaries",
                     expand=c(0,0),
                     breaks=c(0,10,20,30,40,50),
                     limits=c(0,35))+
  scale_fill_brewer(palette="Dark2")+
  scale_alpha_discrete(range=c(0.2,1))+
  guides(fill=guide_legend(title="Repair\nOutcome"))+
  theme_classic(base_size = 18)+ 
  theme(axis.text.x=element_text(size=18, face="bold"),
        axis.text.y=element_text(size=18, face="bold"),
        axis.title.x=element_text(size=18, face="bold"))+
  theme(strip.text.x = element_text(size=20, face="bold"),
        strip.text.y = element_text(size=20, face="bold"))+
  geom_bar(stat="identity", colour="black") +
  scale_x_continuous(limits=c(145,162),
                     breaks=146:161,
                     labels = c("A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C"), 
                     # PolyT1
                     # scale_x_continuous(limits=c(145,162),
                     #                    breaks=146:161,
                     #                    labels = c("C","A","T","T","A"
                     #                               ,"C","C","T","G","T","T","A","T","C","C","C"),
                     # PolyT3 
                     # labels = c("C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C"), 
                     
                     # PolyT4 
                     # labels = c("C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "T", "T", "A", "T", "C", "C", "C"), 
                     
                     expand=c(0,0))
f2



#============================================== Inaccurate Repair and Consistency Plots ==================================================
# csv input - all_SD-MMEJ_consistency, _amplicon_consistency
# csv output - none 
rm(list=ls())
# setwd("/Users/awelterofcollidingmaterials/Desktop/sdmmej-master/Test_Data/PolyT4_output")
library(ggplot2)
all <- read.csv(out_all)
all$PG <- paste("Iw7") 
q<-read.csv("Iw7_amplicon_consistency.csv") # This table is made by hand in excel with all of the consistency stuff for each
# Figure out the SD-MMEJ consistency for each REPAIR_TYPE and put into new table Iw7CR_amplicon_consistency.csv
# Get the percent inaccurate and CONSISTENCY for each REPAIR_TYPE 
# Separate out each individual REPAIR_TYPE
# Look for CONSISTENCY as a function of of percent_inaccurate

CRindel <- subset(all,REPAIR_TYPE=="InDel")
CRmhj <- subset(all,REPAIR_TYPE=="MHJ")
CRabj <- subset(all,REPAIR_TYPE=="ABJ")

CRindelt <- subset(CRindel,CONSISTENCY=="TRUE")
CRindelf <- subset(CRindel,CONSISTENCY=="FALSE")
CRmhjt <- subset(CRmhj,CONSISTENCY=="TRUE")
CRmhjf <- subset(CRmhj,CONSISTENCY=="FALSE")
CRabjt <- subset(CRabj,CONSISTENCY=="TRUE")
CRabjf <- subset(CRabj,CONSISTENCY=="FALSE")

sum(CRindel$percent_inaccurate) #22.47697%
sum(CRindelt$percent_inaccurate) #8.85022%
sum(CRindelf$percent_inaccurate) #13.62675%
sum(CRindelt$percent_inaccurate)/sum(CRindel$percent_inaccurate)*100 #39.37462%

sum(CRmhj$percent_inaccurate) #64.58547%
sum(CRmhjt$percent_inaccurate) #50.5766%
sum(CRmhjf$percent_inaccurate) #14.00887%
sum(CRmhjt$percent_inaccurate)/sum(CRmhj$percent_inaccurate)*100 #78.30956%

sum(CRabj$percent_inaccurate) #11.92767%
sum(CRabjt$percent_inaccurate) #11.30672%
sum(CRabjf$percent_inaccurate) #0.6209485% 
sum(CRabjt$percent_inaccurate)/sum(CRabj$percent_inaccurate)*100 #94.79405%

position <- c("InDel", "MHJ","ABJ")
tbls2 <- aggregate(percent_inaccurate~REPAIR_TYPE+PG,data=all, sum)
con <- subset(all, CONSISTENCY==TRUE) 
tbls3 <- aggregate(percent_inaccurate~PG,data=con, sum)
#correct order
x2 <- subset(tbls2, REPAIR_TYPE=="ABJ")
y2 <- subset(tbls2, REPAIR_TYPE=="MHJ")
z2 <- subset(tbls2, REPAIR_TYPE=="InDel")
xyz2 <- rbind(z2,y2,x2)
levels(xyz2$REPAIR_TYPE)
xyz2$REPAIR_TYPE <- factor(xyz2$REPAIR_TYPE,levels = xyz2$REPAIR_TYPE,ordered = TRUE)

f1 <- ggplot(xyz2, aes(fill=REPAIR_TYPE))+ 
  geom_bar(aes(x=PG, y=percent_inaccurate), position="stack", stat="identity", colour="black") +
  theme_classic(base_size = 9)+
  guides(fill=FALSE)+
  scale_fill_grey(start=0.0, end=1)+
  scale_y_continuous(name = "% Inaccurate Reads",
                     labels=c("0","25","50","75","100"),
                     expand=c(0,0),
                     limits=c(0,120),
                     breaks=c(0,25,50,75,100))+
  scale_x_discrete(name = "Plasmid")+
  theme(axis.text.x=element_text(size=9, face="bold"),
        axis.text.y=element_text(size=9, face="bold"),
        axis.title.x=element_text(size=10, face="bold"),
        axis.title.y=element_text(size=10, face="bold"))
f1
# _InaccurateRepairPlot

f2 <- ggplot(tbls3, aes(fill=PG))+ 
  geom_bar(aes(x=PG, y=percent_inaccurate), width=0.67,stat="identity",colour="black") +
  ylab("% SD-MMEJ consistent")+
  scale_y_continuous(labels=c("0","25","50","75","100"),
                     expand=c(0,0),
                     limits=c(0,101),
                     breaks=c(0,25,50,75,100))+
  xlab("Plasmid")+
  scale_fill_manual(values=c("black", "gray100", "gray75", "gray50"))+
  guides(fill=FALSE)+
  theme_classic(base_size = 9)+
  theme(axis.text.x=element_text(size=9, face="bold"),
        axis.text.y=element_text(size=9, face="bold"),
        axis.title.x=element_text(size=10, face="bold"),
        axis.title.y=element_text(size=10, face="bold"))
f2
# _Consistency

f3<- ggplot(q, aes(fill=PG))+ 
  geom_bar(aes(x=REPAIR_TYPE, y=percent_consistent), width=0.67,position = position_dodge(width=0.75), 
           stat="identity",colour="black") +
  scale_x_discrete(limits=position,
                   expand=c(0,0.5),
                   labels=c("InDel", "MHJ", "ABJ")) +
  ylab("% SD-MMEJ consistent")+
  scale_y_continuous(labels=c("0","25","50","75","100"),
                     expand=c(0,0),
                     limits=c(0,101),
                     breaks=c(0,25,50,75,100))+
  xlab("Repair Junction")+
  guides(fill=guide_legend(title="Plasmid",
                           title.theme = element_text(size=10, face="bold", angle=0),
                           label.theme=element_text(size=8, face="bold", angle=0)))+
  scale_fill_manual(values=c("black", "gray100", "gray75", "gray50"))+
  theme_classic(base_size = 9)+
  theme(legend.position=c(0,1), legend.justification=c(.1,.8),
        legend.key.size=unit(.25, "cm"))+
  theme(axis.text.x=element_text(size=9, face="bold"),
        axis.text.y=element_text(size=9, face="bold"),
        axis.title.x=element_text(size=10, face="bold"),
        axis.title.y=element_text(size=10, face="bold"))
f3
# _Consistency_Breakdown


#======================================== Insertions and Deletion Distro (consistent and inconsistent) ================================================
# csv input - SD-MMEJ_consistency_true_ins_del2
# csv output - none
rm(list=ls())
library(ggplot2)
library(dplyr)
a <- read.csv("Iw7_SD-MMEJ_consistency_true_ins_del2_curated.csv")
# I need to get the mean and median insertion (INSERTION_LENGTH) and deletion (TOTAL_DELETION) lengths 
# Subset out plasmids (above)
ains <- subset(a, CLASS=="Insertion")
adel <- subset(a, CLASS=="deletion")
# Get mean and median for each plasmid 

ains2 <- abs(ains$TOTAL_INSERTION_TRUE)
ains2 <- lapply(ains2, as.numeric)
ains2 <- data.frame(matrix(unlist(ains2), nrow=length(ains2), byrow=T))
ains2 <-cbind(ains, ains2)
colnames(ains2)
ains2 <- rename(ains2, ABS_TOTAL_INSERTION_TRUE = matrix.unlist.ains2...nrow...length.ains2...byrow...T.)
colnames(ains2)[33] <- "ABS_TOTAL_INSERTION_TRUE"

# Now look for SD-MMEJ consistent insertion and deletion lengths
at <- subset(a, CONSISTENCY == "TRUE")
ainst <- subset(at, CLASS=="Insertion")
adelt <- subset(at, CLASS=="deletion")
# Get mean and median for each plasmid 
ains2t <- abs(ainst$TOTAL_INSERTION_TRUE)
ains2t <- lapply(ains2t, as.numeric)
ains2t <- data.frame(matrix(unlist(ains2t), nrow=length(ains2t), byrow=T))
ains2t <-cbind(ainst, ains2t)
colnames(ains2t)
ains2t <- rename(ains2t, ABS_TOTAL_INSERTION_TRUE = matrix.unlist.ains2t...nrow...length.ains2t...byrow...T.)
colnames(ains2t)[33]<-"ABS_TOTAL_INSERTION_TRUE"

ins.mean <- data.frame(PLASMID="Iw7", mean(ains2$ABS_TOTAL_INSERTION_TRUE))
colnames(ins.mean)[2] <- "mean"
ins.median <- data.frame(PLASMID="Iw7", median(ains2$ABS_TOTAL_INSERTION_TRUE))
colnames(ins.median)[2] <- "median"
del.mean <- data.frame(PLASMID="Iw7", mean(adel$TOTAL_DELETION_TRUE))
colnames(del.mean)[2] <- "mean"
del.median <- data.frame(PLASMID="Iw7", median(adel$TOTAL_DELETION_TRUE))
colnames(del.median)[2] <- "median"

# SD-MMEJ Consistent
inst.mean <- data.frame(PLASMID="Iw7", mean(ains2t$ABS_TOTAL_INSERTION_TRUE))
colnames(inst.mean)[2] <- "mean"
inst.median <- data.frame(PLASMID="Iw7", median(ains2t$ABS_TOTAL_INSERTION_TRUE))
colnames(inst.median)[2] <- "median"
delt.mean <- data.frame(PLASMID="Iw7", mean(adelt$TOTAL_DELETION_TRUE))
colnames(delt.mean)[2] <- "mean"
delt.median <- data.frame(PLASMID="Iw7", median(adelt$TOTAL_DELETION_TRUE))
colnames(delt.median)[2] <- "median"

a[, 31][a[, 31] == 0] <- NA
ains2[, 33][ains2[, 33] == 0] <- NA
ains2t[, 33][ains2t[, 33] == 0] <- NA
# First for all junctions
f1 <- ggplot(a)+
  geom_bar(aes(x=TOTAL_DELETION_TRUE, y=percent_inaccurate), position="stack", stat="identity") +
  theme_bw(base_size = 9)+
  facet_wrap(~PG, scales="free_y")+
  geom_vline(aes(xintercept = mean), colour="blue", size=0.75,linetype = "longdash",data=del.mean)+
  geom_vline(aes(xintercept = median), colour="red", size=0.75,linetype = "longdash",data=del.median)+
  scale_y_continuous(name = "Percent Inaccurate Reads",
                     breaks = seq(0,3200, by=10))+  
  scale_x_continuous(name="Deletion Length (bp)",
                     breaks = seq(0, 40, by = 1))+
  theme(axis.text.x=element_text(size=10, face="bold"),
        axis.text.y=element_text(size=10, face="bold"),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold"))+
  theme(strip.text.x = element_text(size=10, face="bold"),
        strip.text.y = element_text(size=10, face="bold"))
(f1)

f2 <- ggplot(ains2)+
  geom_bar(aes(x=ABS_TOTAL_INSERTION_TRUE, y=READS), position="stack", stat="identity") +
  theme_bw(base_size = 9)+
  facet_wrap(~PG, scales="free_y")+
  geom_vline(aes(xintercept = mean), colour="blue", size=0.75,linetype = "longdash",data=ins.mean)+
  geom_vline(aes(xintercept = median), colour="red", size=0.75,linetype = "longdash",data=ins.median)+
  scale_y_continuous(name = "# of Reads",
                     breaks = seq(0,1600, by=2))+  
  scale_x_continuous(name="Insertion Length (bp)",
                     breaks = seq(0,100, by=1))+
  theme(axis.text.x=element_text(size=10, face="bold"),
        axis.text.y=element_text(size=10, face="bold"),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold"))+
  theme(strip.text.x = element_text(size=10, face="bold"),
        strip.text.y = element_text(size=10, face="bold"))
f2

# Now for SD-MMEJ consistent junctions
C1 <- ggplot(adelt)+
  geom_bar(aes(x=TOTAL_DELETION_TRUE, y=percent_inaccurate), position="stack", stat="identity") +
  theme_bw(base_size = 9)+
  facet_wrap(~PG, scales="free_y")+
  geom_vline(aes(xintercept = mean), colour="blue", size=0.75,linetype = "longdash",data=delt.mean)+
  geom_vline(aes(xintercept = median), colour="red", size=0.75,linetype = "longdash",data=delt.median)+
  scale_y_continuous(name = "Percent Inaccurate Reads",
                     breaks = seq(0,3200, by=5))+  
  scale_x_continuous(name="Deletion Length (bp)",
                     breaks = seq(0, 60, by = 1))+
  theme(axis.text.x=element_text(size=10, face="bold"),
        axis.text.y=element_text(size=10, face="bold"),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold"))+
  theme(strip.text.x = element_text(size=10, face="bold"),
        strip.text.y = element_text(size=10, face="bold"))
C1

C2 <- ggplot(ains2t)+
  geom_bar(aes(x=ABS_TOTAL_INSERTION_TRUE, y=percent_inaccurate), position="stack", stat="identity") +
  theme_bw(base_size = 9)+
  facet_wrap(~PG, scales="free_y")+
  geom_vline(aes(xintercept = mean), colour="blue", size=0.75,linetype = "longdash",data=inst.mean)+
  geom_vline(aes(xintercept = median), colour="red", size=0.75,linetype = "longdash",data=inst.median)+
  scale_y_continuous(name = "Percent Inaccurate Reads",
                     breaks = seq(0,1600, by=5))+  
  scale_x_continuous(name="Insertion Length (bp)",
                     breaks = seq(0,60, by=1))+
  theme(axis.text.x=element_text(size=10, face="bold"),
        axis.text.y=element_text(size=10, face="bold"),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold"))+
  theme(strip.text.x = element_text(size=10, face="bold"),
        strip.text.y = element_text(size=10, face="bold"))
C2
# dev.off()

#====================================================== Primer distance and Length =========================================
# csv input - insertion_consistent_analysis
# csv output - 

# rm(list=ls())
m <- read.csv(out_ins_consistency)
m$PG <-paste("Iw7") 
no_trans <- subset(m,is_trans=="no")
trans <- subset(m,is_trans=="trans")
right <- subset(no_trans,side=="RIGHT")
left <- subset(no_trans,side=="LEFT")
right$p1_p2 <- right$p1_start-(right$p2_end+1)
left$p1_p2 <- left$p2_start-(left$p1_end+1)
trans$p1_p2 <- 0
dist <- rbind(right,left,trans)
dist <- dist[order(dist$READS,decreasing = TRUE),]
# no trans
dist.nt <- subset(dist, mechanism!="Trans")

# Need to get the mean and median values for primer length
Iw7a <- subset(dist, PG=="Iw7")
Iw7a <- subset(Iw7a,!(is.na(Iw7a["p1_length"])))

len.mean <- data.frame(PLASMID.x="Iw7", mean(Iw7a$p1_length))
colnames(len.mean)[2] <- "mean"
len.median <- data.frame(PLASMID.x="Iw7",median(Iw7a$p1_length))
colnames(len.median)[2] <- "median"
dist.mean <- data.frame(PLASMID.x="Iw7", mean(Iw7a$p1_p2))
colnames(dist.mean)[2] <- "mean"
dist.median <- data.frame(PLASMID.x="Iw7", median(Iw7a$p1_p2))
colnames(dist.median)[2] <- "median"

dist.nt2 <- abs(dist.nt$p1_p2)
dist.nt2 <- lapply(dist.nt2, as.numeric)
dist.nt2 <- data.frame(matrix(unlist(dist.nt2), nrow=length(dist.nt2), byrow=T))
dist.nt2 <-cbind(dist.nt, dist.nt2)
colnames(dist.nt2)
dist.nt2 <- rename(dist.nt2, p1_p2_abs = matrix.unlist.dist.nt2...nrow...length.dist.nt2...byrow...T.)
colnames(dist.nt2)[77] <- "p1_p2_abs"

f2 <- ggplot(dist.nt2)+
  geom_bar(aes(x=p1_p2_abs, y=percent_inaccurate), position="stack", stat="identity", colour="black", fill="grey50") +
  theme_bw()+
  facet_wrap(~PG, scales="free_y")+
  geom_vline(aes(xintercept = mean), colour="blue", size=0.75,linetype = "longdash",data=dist.mean)+
  geom_vline(aes(xintercept = median), colour="red", size=0.75,linetype = "longdash",data=dist.median)+
  scale_y_continuous(name = "Percent Inaccurate Reads")+  
  scale_x_continuous(name="Distance between Primer Pairs (bp)",
                     breaks = seq(0, 60, by = 1))+
  theme(axis.text.x=element_text(size=10, face="bold"),
        axis.text.y=element_text(size=10, face="bold"),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold"))+
  theme(strip.text.x = element_text(size=10, face="bold"),
        strip.text.y = element_text(size=10, face="bold"))
f2

f1 <- ggplot(dist)+
  geom_bar(aes(x=p2_length, y=percent_inaccurate), position="stack", stat="identity", colour="black", fill="grey50") +
  theme_bw()+
  facet_wrap(~PG, scales="free_y")+
  geom_vline(aes(xintercept = mean), colour="blue", size=0.75,linetype = "longdash",data=len.mean)+
  geom_vline(aes(xintercept = median), colour="red", size=0.75,linetype = "longdash",data=len.median)+
  scale_y_continuous(name = "Percent Inaccurate Reads")+  
  scale_x_continuous(name="Primer Length (bp)",
                     breaks = seq(1, 20, by = 1))+
  theme(axis.text.x=element_text(size=10, face="bold"),
        axis.text.y=element_text(size=10, face="bold"),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold"))+
  theme(strip.text.x = element_text(size=10, face="bold"),
        strip.text.y = element_text(size=10, face="bold"))
f1


# ===================================================== Microhomology length and Use Plots =======================================================
# csv input - SD-MMEJ_consistency_true_ins_del2
# csv output - none 
rm(list=ls())
library(grid)
library(ggplot2)
# setwd("/Users/awelterofcollidingmaterials/Desktop/sdmmej-master/Test_Data/FlexLoopSBMH/FlexLoopSBMH_output")
r<-read.csv(out_consistency_true_ins)
# setwd("/Users/awelterofcollidingmaterials/Desktop/Iw7R_HiSeq/Flexible Loop Analysis/")
# f<-read.csv("flexloops.csv")
# f$percent_inaccurate <- f$READS/sum(f$READS)*100
# sum(f$percent_inaccurate)
# fL<-subset(f, side=="LEFT")
# fR<-subset(f, side=="RIGHT")

# Investigating microhomologies produced during insertion events. I should subset out L and R repeat motifs and then do the MH analysis
r[, 22][r[, 22] == 0] <- NA
T1 <- ggplot(r)+
  geom_bar(aes(x=MH_Length, y=percent_inaccurate), position="stack", stat="identity", colour="black", fill="grey50") +
  theme_bw()+
  facet_wrap(~PG, scales="free_y")+
  # geom_vline(aes(xintercept = mean), colour="blue", size=0.75,linetype = "longdash",data=len.mean)+
  # geom_vline(aes(xintercept = median), colour="red", size=0.75,linetype = "longdash",data=len.median)+
  scale_y_continuous(name = "Percent Inaccurate Reads")+
  scale_x_continuous(name="Microhomology (bp)",
                     breaks = seq(1, 100, by = 1))+
  theme(axis.text.x=element_text(size=10, face="bold"),
        axis.text.y=element_text(size=10, face="bold"),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold"))+
  theme(strip.text.x = element_text(size=10, face="bold"),
        strip.text.y = element_text(size=10, face="bold"))
T1

# (hsb7 <- hsb2.small[hsb2.small$id %in% c(12, 48, 86, 11, 20, 195), ])


a<-r[r$MICROHOMOLOGY %in% c("A","G","C","T","CC","CCT","TA","TC"),] #you'll want to use this sometimes
T1 <- ggplot(a)+
  geom_histogram(aes(x=MICROHOMOLOGY, y=percent_inaccurate), position="stack", stat="identity", colour="black", fill="grey50") +
  theme_bw()+
  # facet_wrap(~PG, scales="free_y")+
  # geom_vline(aes(xintercept = mean), colour="blue", size=0.75,linetype = "longdash",data=len.mean)+
  # geom_vline(aes(xintercept = median), colour="red", size=0.75,linetype = "longdash",data=len.median)+
  # ggtitle("Left Synthesis")+
  scale_y_continuous(name = "Percent Inaccurate")+  
  scale_x_discrete(name="Microhomology")+
  theme(plot.title = element_text(color="grey28", size=12),
        axis.text.x=element_text(size=10, face="bold"),
        axis.text.y=element_text(size=10, face="bold"),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold"))+
  theme(strip.text.x = element_text(size=10, face="bold"),
        strip.text.y = element_text(size=10, face="bold"))
T1


a<-read.csv(out_break)
# MH Plot (Deletions)
# a2<-subset(r, MICROHOMOLOGY!=0)
a$percent_inaccurate <- a$READS/sum(a$READS)*100
sum(a$percent_inaccurate)
aL<-subset(a, side=="LEFT")
aR<-subset(a, side=="RIGHT")
# Microhomology Plot (insertions)
T1 <- ggplot(aL)+
  geom_histogram(aes(x=insMH, y=percent_inaccurate), position="stack", stat="identity", colour="black", fill="grey50") +
  theme_bw()+
  # facet_wrap(~PG, scales="free_y")+
  # geom_vline(aes(xintercept = mean), colour="blue", size=0.75,linetype = "longdash",data=len.mean)+
  # geom_vline(aes(xintercept = median), colour="red", size=0.75,linetype = "longdash",data=len.median)+
  ggtitle("Left Synthesis")+
  scale_y_continuous(name = "Percent Inaccurate")+  
  scale_x_discrete(name="Microhomology")+
  theme(plot.title = element_text(color="grey28", size=12),
        axis.text.x=element_text(size=10, face="bold"),
        axis.text.y=element_text(size=10, face="bold"),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold"))+
  theme(strip.text.x = element_text(size=10, face="bold"),
        strip.text.y = element_text(size=10, face="bold"))
T1

T2 <- ggplot(aR)+
  geom_histogram(aes(x=insMH, y=percent_inaccurate), position="stack", stat="identity", colour="black", fill="grey50") +
  theme_bw()+
  # facet_wrap(~PG, scales="free_y")+
  # geom_vline(aes(xintercept = mean), colour="blue", size=0.75,linetype = "longdash",data=len.mean)+
  # geom_vline(aes(xintercept = median), colour="red", size=0.75,linetype = "longdash",data=len.median)+
  ggtitle("Right Synthesis")+
  scale_y_continuous(name = "Percent Inaccurate")+  
  scale_x_discrete(name="Microhomology")+
  theme(plot.title = element_text(color="grey28", size=12),
        axis.text.x=element_text(size=10, face="bold"),
        axis.text.y=element_text(size=10, face="bold"),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold"))+
  theme(strip.text.x = element_text(size=10, face="bold"),
        strip.text.y = element_text(size=10, face="bold"))

T2
vp1 <- viewport(layout=grid.layout(2,1), width = 1, height = 1/2,y=.75)
vp2 <- viewport(layout=grid.layout(2,1),width = 1, height = 1/2,y=0.25)
vp1
print(T1, vp = vp1)
print(T2, vp = vp2)

# FREE ENERGY DISTRIBUTION
# FE <-read.csv("PolyGTFreeEnergy.csv")
# 
# T1<-boxplot(FE$PolyG,FE$PolyT,ylim = c(0, -18), 
#             main="Secondary Structure Stability",
#             names=c("Poly G Loop", "Poly T Loop"),
#             xlab="Plasmid",
#             ylab="DeltaG",
#             col="#69b3a2",
#             border="brown"
# )
# T1$stats
# # p-value = 0.013 
# 
# t.test(FE$PolyG,FE$PolyT)
# 
# 
# nbGroup <- nlevels(FE$PolyT)
# text( 
#   x=c(1:nbGroup), 
#   y=T1$stats[nrow(T1$stats),] + 0.5, 
#   paste("n = ",table(FE$PolyT),sep="")  
# )