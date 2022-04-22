#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressWarnings(suppressMessages({
  library(seqinr)
  library(ggplot2)
  library(grid)
  library(gtable)
  library(stringr)
  library(dplyr)
  library(Biostrings)
  library(tidyverse)
  library(reshape)
  library(plyr)
}))


##test if there is at least 2 argument: if not, return an error
if (length(args)<2) {
  stop("Usage: Rscript SDMMEJ_Plots.R outdir plasmid", call.=FALSE)
} else if (length(args)== 2) {
  outdir=args[1]
  plasmid=args[2]
}

#In order to debug, uncomment this and comment the command line information
outdir="/Users/rbator01/Library/CloudStorage/Box-Box/bioinformatics_research_technology/rt_bioinformatics_consultations/mcvey_lab_rt_bioinformatics/terrence_dna_repair/HiSeq_CRISPR_Data/renamed_files/R1_output/"
plasmid="R1"

#-----------------------Read In Data Frames--------------------------------------

## Some of these files may be missing from the pipeline, in that case, we need to create an empty dataframe so the code will run
reclassData = read.csv(paste0(outdir, "/", plasmid, "_reclassified.csv"), row.names = 1)

## Insertion Files
insertion_file=paste0(outdir, "/", plasmid,"_insertion_insertion_consistency2.csv")
ins_consistent_colnames=c("ID","DR_START","DR_END","RC_START","RC_END","consistency","RECONSTRUCTED_SEQ",
                         "left_del","right_del","del_seq","insertion","plasmid","DRmotif_length",
                         "RCmotif_length","Loop-out","Snap-back")

if (file.exists(insertion_file)){
  insertionsConsistent = read.csv(insertion_file, row.names=1)
}else{
  insertionsConsistent = data.frame(matrix(nrow=1, ncol=length(ins_consistent_colnames)))
  colnames(insertionsConsistent) = ins_consistent_colnames
}
complex_file=paste0(outdir, "/", plasmid, "_complex_insertion_consistency2.csv")

if (file.exists(complex_file)){
  complexConsistent = read.csv(complex_file, row.names = 1)
  
}else{
  complexConsistent = data.frame(matrix(nrow=1, ncol=length(ins_consistent_colnames)))
  colnames(complexConsistent) = ins_consistent_colnames
  
}


## Deletion Files
del_file_1=paste0(outdir, "/", plasmid, "_deletion_consistency_log_subset.csv")
del_file_2=paste0(outdir, "/", plasmid, "_deletion_consistency_table.txt")

if(file.exists(del_file_1)){
  del = read.csv(del_file_1)
}else{
  colnames=c("Sample ID","Deletion Length","Repair Type","Mechanism","Motif to Break","Motif to Deletion","P1 to Break","P1 to Deletion","P2 to Break","P2 to Deletion","P1 to P2","Motif Length","Break Side",
             "Deletion to MH","Motif Sequence")
  del=data.frame(matrix(nrow=1, ncol=length(colnames)))
  colnames(del) = colnames
}

if(file.exists(del_file_2)){
  DeletionData = read.csv(del_file_2, header=T, sep="\t")
}else{
  colnames=c("PLASMID","ID","GENOTYPE","REPAIR_TYPE","SEQ","TOTAL_DELETION","LEFT_DEL_INDEX","RIGHT_DEL_INDEX","CONSISTENCY")
  DeletionData = data.frame(matrix(nrow=1, ncol=length(colnames)))
  colnames(DeletionData) = colnames
}


#----------------------Basic Information Input-----------------------------------
names(reclassData)[names(reclassData) == "ALIGNED_SEQ"] <- "RECONSTRUCTED_SEQ"

#### get the length of the final sequence
referenceSeq = subset(reclassData, CLASS=="exact")
referencetxt = referenceSeq$RECONSTRUCTED_SEQ
referenceSplit = unlist(strsplit(as.character(referencetxt), split=""))
lengthRef = nchar(as.character(referencetxt))

BreakPointFromLeft =  referenceSeq$DISTANCE_FROM_BREAK_RIGHT #Input Break Point From Left used in Hi-FiBR
BreakPointFromRight =  referenceSeq$DISTANCE_FROM_BREAK_LEFT #Input Break Point From Right used in Hi-FiBR

#-----------------------------Data Manipulation for further SD-MMEJ Analysis-------------------------------------------
CombinedCurated = subset(reclassData, READS>9) #eliminates all lines with less reads than what we choose to specify 
CombinedCurated$percent = CombinedCurated$READS/sum(CombinedCurated$READS)*100
inaccurate = subset(CombinedCurated, CLASS!="exact") #removes the "exact" sequence to narrow down to inaccurate repair events only
accurate = subset(CombinedCurated,CLASS=="exact")
accurate$percent_inaccurate <- "NA"
inaccurate$percent_inaccurate <- inaccurate$READS/sum(inaccurate$READS)*100
CombinedCurated = rbind(accurate,inaccurate)

#reorganizing data by repair event classification
deletions = subset(CombinedCurated, CLASS_final=="deletion")
head(deletions)
complex = subset(CombinedCurated, CLASS_final=="complex")
insertions = subset(CombinedCurated, CLASS_final=="insertion")
exact = subset(CombinedCurated, CLASS_final=="exact")
CombinedCurated=rbind(deletions, complex, insertions, exact)

write.csv(CombinedCurated, paste0(outdir, "/", "table_outputs/",plasmid,"_combined_curated.csv"))

complexInsertion=rbind(complex, insertions)

InsComConsistent = rbind(insertionsConsistent, complexConsistent)
InsComConsistent$RIGHT_DEL <- lengthRef -(nchar(as.character(InsComConsistent$RECONSTRUCTED_SEQ))-InsComConsistent$right_del+1) 
InsComMerge = merge(InsComConsistent, complexInsertion, "RECONSTRUCTED_SEQ", all=FALSE)

# these tables can be empty
#try( {del$SEQ = toupper(del$SEQ)}, silent=TRUE)
del$SEQ = toupper(del$SEQ)
names(del)[names(del)=="SEQ"] = "RECONSTRUCTED_SEQ"

deletionsMerged = merge(del, deletions, "RECONSTRUCTED_SEQ", all=FALSE )

write.csv(deletionsMerged, paste(outdir, "/", "table_outputs/",plasmid,"_deletions_table.csv", sep=""))

InsertsComplex = InsComMerge[,18:32]
InsertsComplex$RECONSTRUCTED_SEQ = InsComMerge$RECONSTRUCTED_SEQ
InsertsComplex$READS = InsComMerge$READS
InsertsComplex$MICROHOMOLOGY = InsComMerge$MICROHOMOLOGY 
InsertsComplex$MH_Length = InsComMerge$MH_Length
InsertsComplex$NUMBER_OF_ALIGNMENTS  = InsComMerge$NUMBER_OF_ALIGNMENTS
InsertsComplex$percent = InsComMerge$percent
InsertsComplex$percent_inaccurate = InsComMerge$percent_inaccurate

# these tables can be empty
InsertsComplex$REPAIR_TYPE = "InDel"
InsertsComplex$CLASS = "Insertion"

# try( {InsertsComplex$REPAIR_TYPE = "InDel"}, silent=TRUE)
# try( {InsertsComplex$CLASS = "Insertion"}, silent=TRUE)

InsertsComplex$left_del = InsComMerge$left_del 
InsertsComplex$right_del = InsComMerge$right_del
InsertsComplex$CONSISTENCY = InsComMerge$consistency

del1 <- cbind(deletionsMerged[,10:24])

del1$INSERTED_SEQ <- NA
del1$RECONSTRUCTED_SEQ <- deletionsMerged$RECONSTRUCTED_SEQ
del1$READS <- deletionsMerged$READS
del1$MICROHOMOLOGY <- deletionsMerged$MICROHOMOLOGY 
del1$MH_Length <- deletionsMerged$MH_Length
del1$NUMBER_OF_ALIGNMENTS  <- deletionsMerged$NUMBER_OF_ALIGNMENTS
del1$percent <- deletionsMerged$percent
del1$percent_inaccurate <- deletionsMerged$percent_inaccurate
del1$REPAIR_TYPE <- deletionsMerged$REPAIR_TYPE
del1$CLASS = deletionsMerged$CLASS_final
deletions <- cbind(del1, deletionsMerged[,7:9])
names(deletions)
names(InsertsComplex)

## can't do this if one is empty, there aren't the same colnames
names(deletions) <- names(InsertsComplex)


all <- rbind(InsertsComplex,deletions)
view(all)

all$plasmid = plasmid

write.csv(all,paste(outdir, "/", "table_outputs/",plasmid,"_all_SD-MMEJ_consistency.csv", sep=""))

#Calculating efficiency of inaccurate repair events returned
ExactReads = exact$READS
InaccurateReads = sum(all$READS)
AllReads = sum(ExactReads, InaccurateReads)
PercentInaccurateReads = InaccurateReads/AllReads*100
PercentInaccurateReadsTable = as.data.frame(PercentInaccurateReads)
PercentInaccurateReadsTable$plasmid=plasmid
write.csv(PercentInaccurateReads, paste(outdir, "/", "table_outputs/", plasmid, "_percent_inaccurate_reads.csv", sep=""))

PercentInaccurateReadsPlot = ggplot(PercentInaccurateReadsTable, aes(fill=plasmid))+geom_bar(aes(x=plasmid ,y=PercentInaccurateReads), width=0.67, stat="identity", colour="black")+
  ylab("% Inaccurate Repair Reads")+scale_y_continuous(expand=c(0,0),limits=c(0,PercentInaccurateReadsTable$PercentInaccurateReads+5))+xlab("plasmid")+
  scale_fill_manual(values=c("black", "gray100", "gray75", "gray50"))+guides(fill="none")+theme_classic(base_size = 9)+
  theme(axis.text.x=element_text(size=9, face="bold"),axis.text.y=element_text(size=9, face="bold"),axis.title.x=element_text(size=10, face="bold"),
        axis.title.y=element_text(size=10, face="bold"))
PercentInaccurateReadsPlot
pdf(paste(outdir, "/", "plots/", plasmid, "_Percent_Inaccurate_Repair_Reads_Plot.pdf", sep=""), width=5)
PercentInaccurateReadsPlot
dev.off()

# Subset combined curated file by consistent complex and insertions for repeat motif and primer plots
insertions2 = subset(CombinedCurated, CLASS_final=="insertion")
complex2 = subset(CombinedCurated, CLASS_final=="complex")
CombCurInsertions = rbind(insertions2,complex2)
InsComConsistentTrue = subset(InsComConsistent, consistency==TRUE)
ConsistentInsertions = merge(CombCurInsertions, InsComConsistentTrue, by="RECONSTRUCTED_SEQ") #combines all the consistent insertions by SEQ so that READS are now attributed to other details

#----------------Direct Repeat Plot Data Manipulation-----------------------------------------
ConsistentInsertions = ConsistentInsertions[!is.na(ConsistentInsertions$DR_START),]
ConsistentInsertions$SIDE = ifelse(ConsistentInsertions$DR_END > ConsistentInsertions$right_del, "RIGHT", "LEFT")
LoopOut = as.data.frame(ConsistentInsertions$Loop.out)
LoopOut = as.data.frame(sapply(LoopOut, gsub, pattern="-", replacement=""))

# can be empty
try({colnames(LoopOut)[1]="RM"}, silent=TRUE)

ConsistentInsertions = cbind(ConsistentInsertions, LoopOut)
ConsistentInsertions = ConsistentInsertions[!is.na(ConsistentInsertions$RM),]
ConsistentInsertions$insertion_length = nchar(as.character(ConsistentInsertions$insertion))

ConsistentInsertions$DR_START_TRUE = ifelse(ConsistentInsertions$SIDE == "RIGHT", 
         (lengthRef - ConsistentInsertions$READ_LENGTH + ConsistentInsertions$DR_START + ConsistentInsertions$insertion_length), 
         ConsistentInsertions$DR_START)
# if right side, then, (Reference sequence length - Read Length + DR_START + insertion length)), 
#     if left side it takes the DR_START value. Track the parenthesis if the values are returning negative

ConsistentInsertions$DR_END_TRUE = ifelse(ConsistentInsertions$SIDE=="RIGHT",
         (lengthRef - ConsistentInsertions$READ_LENGTH + ConsistentInsertions$DR_END + ConsistentInsertions$insertion_length),
         ConsistentInsertions$DR_END)

#if right side, then, (Reference sequence length - Read Length + DR_END + insertion length)), 
#     if left side, it take the DR_END value. Track the parenthesis if the values are returning negative
            #What these lines just above are doing is correcting for previous scripts. They are taking into account that the    #
            #repeat motif that exists on the right of the break is actually the second instance of the repeated sequence.       #
            #This is subtracting the read length from the reference, giving a negative number, then adding back where previous  #
            #scripts thought the Direct Repeat Motifs stat/end are, plus the insertion length.                                  #

ConsistentInsertions$motif_ID = paste(ConsistentInsertions$ID, ConsistentInsertions$DR_START_TRUE, ConsistentInsertions$DR_END_TRUE, sep="-")
ConsistentInsertions = ConsistentInsertions[!duplicated(ConsistentInsertions$motif_ID),] #Removes duplicate repeat motif IDs
            #The duplicates that get removed are due to SD-MMEJ events that could occur by either mechanism, therefore they     #
            #are included in both the direct repeat motifs and the reverse complement repeat motifs.                            #

ConsistentInsertionsDRR = subset(ConsistentInsertions, SIDE == "RIGHT")
ConsistentInsertionsDRL = subset(ConsistentInsertions, SIDE == "LEFT")
RPrimerPlotPoints = c() # Creates new data frame "RPrimerPlotPoints" for the for loop below
LPrimerPlotPoints = c() # Creates new data frame "LPrimerPlotPoints" for the for loop below

#----------------Direct Repeat- Right Side-----------------------------------------
if (nrow(ConsistentInsertionsDRR) > 0){
  for (i in 1:nrow(ConsistentInsertionsDRR)){ 
    Ins <- as.character(ConsistentInsertionsDRR[i,35])    # creates insertion sequence region to search using "insertion" column
    dr1 <- as.data.frame(str_locate_all(ConsistentInsertionsDRR[i,43], Ins)) # search deletion sequence for string using "RM" column
    if (nrow(dr1)!=1){
      dr2 <- dr1
      while (nrow(dr2)!=1){
        for (j in 1:nrow(dr1)){
          if (dr1[j,2]==nchar(allowNA=FALSE, as.character(ConsistentInsertionsDRR[i,43]))){ #RM column
            break
          }else{
            prim <- paste(substring(ConsistentInsertionsDRR[i,43],1,dr1[j,1]-1), substring(ConsistentInsertionsDRR[i,43],dr1[j,2]+1,nchar(as.character(ConsistentInsertionsDRR[i,43]))),sep="")
            s <- dr1[j,1]-1
            der <- substring(as.character(ConsistentInsertionsDRR[i,34]), first=ConsistentInsertionsDRR[i,32]-(s-1), last=ConsistentInsertionsDRR[i,33]+(nchar(as.character(ConsistentInsertionsDRR[i,43]))-dr1[j,2]))
            #delSeq column                              #left_del                         #right_del
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
      RPrimerPlotPoints1 <- cbind(ConsistentInsertionsDRR[i,47], dr1[k,]) # once we know which one we want, we put that dr1 row
      #motif_ID
      RPrimerPlotPoints <- rbind(RPrimerPlotPoints,RPrimerPlotPoints1)
    } else {
      RPrimerPlotPoints1 <- cbind(ConsistentInsertionsDRR[i,47], dr1)
      RPrimerPlotPoints <- rbind(RPrimerPlotPoints,RPrimerPlotPoints1)
    }
  }
  colnames(RPrimerPlotPoints)[1] = "motif_ID"
  RPrimerPlotPoints = RPrimerPlotPoints[!duplicated(RPrimerPlotPoints$motif_ID),] #Removes duplicate repeat motif IDs
  #This shouldn't remove anything because they were removed above, but just making sure so nothing is double counted    #
  
  RPrimerPlotPoints$is_trans = ifelse(ConsistentInsertionsDRR$DR_START_TRUE<BreakPointFromLeft,
                                      ifelse(ConsistentInsertionsDRR$DR_START_TRUE>BreakPointFromLeft, "trans", "no"), "no")
  #I don't understand this ifelse statement above, its illogical and will always lead to no       #
  
  RPrimerPlotPoints$p1_start = ifelse(RPrimerPlotPoints$is_trans=="trans",
                                      ConsistentInsertionsDRR$DR_START_TRUE, (ConsistentInsertionsDRR$DR_END_TRUE+1)-
                                        (nchar(as.character(ConsistentInsertionsDRR$RM))-RPrimerPlotPoints$end))
  RPrimerPlotPoints$p1_end = ifelse(RPrimerPlotPoints$is_trans=="trans",
                                    ConsistentInsertionsDRR$DR_START_TRUE+(RPrimerPlotPoints$start-2),
                                    ConsistentInsertionsDRR$DR_END_TRUE)
  RPrimerPlotPoints$p1_length = 1+(RPrimerPlotPoints$p1_end - RPrimerPlotPoints$p1_start)
  RPrimerPlotPoints$mh_start = ifelse(RPrimerPlotPoints$is_trans=="trans", 
                                      (ConsistentInsertionsDRR$DR_END_TRUE+1)-(nchar(as.character(ConsistentInsertionsDRR$RM))-RPrimerPlotPoints$end), 
                                      ConsistentInsertionsDRR$DR_START_TRUE)
  RPrimerPlotPoints$mh_end = ifelse(RPrimerPlotPoints$is_trans=="trans",
                                    ConsistentInsertionsDRR$DR_END_TRUE, ConsistentInsertionsDRR$DR_START_TRUE+(RPrimerPlotPoints$start-2))
  RPrimerPlotPoints$mh_length = 1+(RPrimerPlotPoints$mh_end-RPrimerPlotPoints$mh_start)
  RPrimerPlotPoints$ins_start = ConsistentInsertionsDRR$DR_START_TRUE+(RPrimerPlotPoints$start-1)
  RPrimerPlotPoints$ins_end = (ConsistentInsertionsDRR$DR_END_TRUE)-(nchar(as.character(ConsistentInsertionsDRR$RM))-RPrimerPlotPoints$end)
  RPrimerPlotPoints$ins_length = 1+(RPrimerPlotPoints$ins_end-RPrimerPlotPoints$ins_start)
  RPrimerPlotPoints$p2_start = ifelse(RPrimerPlotPoints$is_trans=="trans",
                                      RPrimerPlotPoints$p1_start, ConsistentInsertionsDRR$RIGHT_DEL+1)
  RPrimerPlotPoints$p2_end = ifelse(RPrimerPlotPoints$is_trans=="trans",
                                    RPrimerPlotPoints$p1_end, (ConsistentInsertionsDRR$RIGHT_DEL+1)+(RPrimerPlotPoints$p1_length-1))
  RPrimerPlotPoints$p2_length = 1+(RPrimerPlotPoints$p2_end-RPrimerPlotPoints$p2_start)
  ConsistentInsertionsDRR = merge(ConsistentInsertionsDRR, RPrimerPlotPoints, by="motif_ID")
  ConsistentInsertionsDRR$p_ID = paste(ConsistentInsertionsDRR$p1_start, ConsistentInsertionsDRR$p1_end, ConsistentInsertionsDRR$p2_start, ConsistentInsertionsDRR$p2_end, sep="-")
  DRRagg = aggregate(READS~p_ID, data=ConsistentInsertionsDRR, sum)
  DRRtable = as.data.frame(table(ConsistentInsertionsDRR$p_ID))
  DRRagg = merge(DRRagg, DRRtable, by.x="p_ID", by.y="Var1")
  
  ConsistentInsertionsDRR = merge(DRRagg, ConsistentInsertionsDRR, by="p_ID")
  ConsistentInsertionsDRRsimplified = ConsistentInsertionsDRR[!duplicated(ConsistentInsertionsDRR$p_ID),]
  ConsistentInsertionsDRRsimplified$p2_start2 = ConsistentInsertionsDRRsimplified$p2_start-0.5
  ConsistentInsertionsDRRsimplified$p2_end2 = ConsistentInsertionsDRRsimplified$p2_end+0.5
  ConsistentInsertionsDRRsimplified = ConsistentInsertionsDRRsimplified[, -c(5:27), drop=FALSE]
  colnames(ConsistentInsertionsDRRsimplified)[2] = "READS"
  }else {
    colnames=c('p_ID','READS','Freq.x','Freq.y','NUMBER_OF_ALIGNMENTS','MISMATCH_PERCENTAGE_TO_RECONSTRUCTED','CLASS_final',
               'percent','percent_inaccurate','ID','DR_START','DR_END','RC_START','RC_END','consistency','left_del','right_del','del_seq','insertion','plasmid',
               'DRmotif_length','RCmotif_length','Loop.out','Snap.back','RIGHT_DEL','SIDE','RM','insertion_length','DR_START_TRUE','DR_END_TRUE','start.x','end.x',
               'is_trans.x','p1_start.x','p1_end.x','p1_length.x','mh_start.x','mh_end.x','mh_length.x','ins_start.x','ins_end.x','ins_length.x','p2_start.x','p2_end.x',
               'p2_length.x','start.y','end.y','is_trans.y','p1_start.y',
               'p1_end.y','p1_length.y','mh_start.y','mh_end.y','mh_length.y','ins_start.y','ins_end.y','ins_length.y','p2_start.y','p2_end.y','p2_length.y')
    ConsistentInsertionsDRRsimplified = data.frame(matrix(nrow=0, ncol=length(colnames)))
    colnames(ConsistentInsertionsDRRsimplified) = colnames
  }                         
              #Deletes some unnecessary columns to make the data frame easier to look through. Changes READS.x to READS   #
              #There's more that can be removed still but leaving it at this point for now.                                #

#----------------Direct Repeat - Left Side-----------------------------------------
if (nrow(ConsistentInsertionsDRL) > 0){
  for (i in 1:nrow(ConsistentInsertionsDRL)){
    Ins <- as.character(ConsistentInsertionsDRL[i,35])    # create insertion sequence region to search
    #insertion
    dr1 <- as.data.frame(str_locate_all(ConsistentInsertionsDRL[i,43], Ins)) # search deletion sequence for string #43 = RM
    if (nrow(dr1)!=1){
      dr2 <- dr1
      while (nrow(dr2)!=1){
        for (j in 1:nrow(dr1)){
          if (dr1[j,2]==nchar(allowNA = TRUE,as.character(ConsistentInsertionsDRL[i,43]))){
            break
          }else{
            prim <- paste(substring(ConsistentInsertionsDRL[i,43],1,dr1[j,1]-1), substring(ConsistentInsertionsDRL[i,43],dr1[j,2]+1,nchar(as.character(ConsistentInsertionsDRL[i,43]))),sep="")
            s <- dr1[j,1]-1
            der <- substring(as.character(ConsistentInsertionsDRL[i,34]), first=ConsistentInsertionsDRL[i,32]-(s-1), last=ConsistentInsertionsDRL[i,33]+(nchar(as.character(ConsistentInsertionsDRL[i,43]))-dr1[j,2]))
            #del_seq                              left_del                               right_del
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
      LPrimerPlotPoints1 <- cbind(ConsistentInsertionsDRL[i,47], dr1[k,]) # once we know which one we want, we put that dr1 row
      #motif_ID
      LPrimerPlotPoints <- rbind(LPrimerPlotPoints,LPrimerPlotPoints1)
    } else {
      LPrimerPlotPoints1 <- cbind(ConsistentInsertionsDRL[i,47], dr1)
      #motif_ID
      LPrimerPlotPoints <- rbind(LPrimerPlotPoints,LPrimerPlotPoints1)
    }
  }
  
  colnames(LPrimerPlotPoints)[1] = "motif_ID"
  LPrimerPlotPoints$is_trans = ifelse(ConsistentInsertionsDRL$DR_END_TRUE<BreakPointFromLeft,
                                      ifelse(ConsistentInsertionsDRL$DR_END_TRUE>BreakPointFromLeft, "trans", "no"), "no")
  #I don't understand this ifelse statement above, its illogical and will always lead to no       #
  
  LPrimerPlotPoints$p1_start = ifelse(LPrimerPlotPoints$is_trans=="trans",
                                      (ConsistentInsertionsDRL$DR_END_TRUE+1)-
                                        (nchar(as.character(ConsistentInsertionsDRR$RM))-LPrimerPlotPoints$end), ConsistentInsertionsDRL$DR_START_TRUE)
  LPrimerPlotPoints$p1_end = ifelse(LPrimerPlotPoints$is_trans=="trans",
                                    ConsistentInsertionsDRL$DR_END_TRUE, ConsistentInsertionsDRL$DR_START_TRUE+(LPrimerPlotPoints$start-2))
  LPrimerPlotPoints$p1_length = 1+(LPrimerPlotPoints$p1_end - LPrimerPlotPoints$p1_start)
  LPrimerPlotPoints$mh_start = ifelse(LPrimerPlotPoints$is_trans=="trans", 
                                      ConsistentInsertionsDRL$DR_START_TRUE, 
                                      (ConsistentInsertionsDRL$DR_END_TRUE+1)-(nchar(as.character(ConsistentInsertionsDRL$RM))-LPrimerPlotPoints[,3]))
  LPrimerPlotPoints$mh_end = ifelse(LPrimerPlotPoints$is_trans=="trans",
                                    ConsistentInsertionsDRL$DR_START_TRUE+(LPrimerPlotPoints$start-2),
                                    ConsistentInsertionsDRL$DR_END_TRUE)
  LPrimerPlotPoints$mh_length = 1+(LPrimerPlotPoints$mh_end-LPrimerPlotPoints$mh_start)
  LPrimerPlotPoints$ins_start = ConsistentInsertionsDRL$DR_START_TRUE+(LPrimerPlotPoints$start-1)
  LPrimerPlotPoints$ins_end = (ConsistentInsertionsDRL$DR_END_TRUE)-(nchar(as.character(ConsistentInsertionsDRL$RM))-LPrimerPlotPoints$end)
  LPrimerPlotPoints$ins_length = 1+(LPrimerPlotPoints$ins_end-LPrimerPlotPoints$ins_start)
  LPrimerPlotPoints$p2_start = ifelse(LPrimerPlotPoints$is_trans=="trans",
                                      LPrimerPlotPoints$p1_start, (ConsistentInsertionsDRL$left_del)-(LPrimerPlotPoints$p1_length-1))
  LPrimerPlotPoints$p2_end = ifelse(LPrimerPlotPoints$is_trans=="trans",
                                    LPrimerPlotPoints$p1_end, ConsistentInsertionsDRL$left_del)
  LPrimerPlotPoints$p2_length = 1+(LPrimerPlotPoints$p2_end-LPrimerPlotPoints$p2_start)
  ConsistentInsertionsDRL = merge(ConsistentInsertionsDRL, LPrimerPlotPoints, by="motif_ID")
  ConsistentInsertionsDRL$p_ID = paste(ConsistentInsertionsDRL$p1_start, ConsistentInsertionsDRL$p1_end, ConsistentInsertionsDRL$p2_start, ConsistentInsertionsDRL$p2_end, sep="-")
  DRLagg = aggregate(READS~p_ID, data=ConsistentInsertionsDRL, sum)
  DRLtable = as.data.frame(table(ConsistentInsertionsDRL$p_ID))
  DRLagg = merge(DRLagg, DRLtable, by.x="p_ID", by.y="Var1")
  ConsistentInsertionsDRL = merge(DRLagg, ConsistentInsertionsDRL, by="p_ID")
  ConsistentInsertionsDRLsimplified = ConsistentInsertionsDRL[!duplicated(ConsistentInsertionsDRL$p_ID),] #removes duplicated p2_IDs for use in repeat motif plots
  ConsistentInsertionsDRLsimplified$p2_start2 = ConsistentInsertionsDRLsimplified$p2_start-0.5
  ConsistentInsertionsDRLsimplified$p2_end2 = ConsistentInsertionsDRLsimplified$p2_end+0.5
  ConsistentInsertionsDRLsimplified = ConsistentInsertionsDRLsimplified[, -c(5:27), drop=FALSE]
  colnames(ConsistentInsertionsDRLsimplified)[2] = "READS"
}else{
  colnames=colnames(ConsistentInsertionsDRRsimplified)
  ConsistentInsertionsDRLsimplified = data.frame(matrix(nrow=0, ncol=length(colnames)))
  colnames(ConsistentInsertionsDRLsimplified) = colnames
}
              #Deletes some unnecessary columns to make the data frame easier to look through. Changes READS.x to READS   #
              #Theres more that can be removed still but leaving it at this point for now.                                #

#----------------Reverse Complement Plot Data Manipulation-----------------------------------------
ConsistentInsertionsRC = merge(CombCurInsertions, InsComConsistentTrue, by="RECONSTRUCTED_SEQ")
ConsistentInsertionsRC = ConsistentInsertionsRC[!is.na(ConsistentInsertionsRC$RC_START),]
ConsistentInsertionsRC$SIDE = ifelse(ConsistentInsertionsRC$RC_END>ConsistentInsertionsRC$right_del, "RIGHT", "LEFT")
SnapBack = as.data.frame(ConsistentInsertionsRC$Snap.back) #makes dataframe of just the snapback column
SnapBack = as.data.frame(sapply(SnapBack,gsub,pattern="-", replacement=""))

# can be empty
try({colnames(SnapBack)[1]="RM"}, silent=TRUE)

ConsistentInsertionsRC = cbind(ConsistentInsertionsRC, SnapBack)
ConsistentInsertionsRC$insertion_length = nchar(as.character(ConsistentInsertionsRC$insertion))
ConsistentInsertionsRC$RC_START_TRUE = ifelse(ConsistentInsertionsRC$SIDE=="RIGHT",
                                              (lengthRef-ConsistentInsertionsRC$READ_LENGTH + ConsistentInsertionsRC$RC_START + ConsistentInsertionsRC$insertion_length),
                                              ConsistentInsertionsRC$RC_START)
# if right side, then, (Reference sequence length - Read Length + RC_START + insertion length), 
#     if left side it takes the RC_START value. Track the parenthesis if the values are returning negative
ConsistentInsertionsRC$RC_END_TRUE = ifelse(ConsistentInsertionsRC$SIDE=="RIGHT",
                                            (lengthRef-ConsistentInsertionsRC$READ_LENGTH + ConsistentInsertionsRC$RC_END + ConsistentInsertionsRC$insertion_length),
                                            ConsistentInsertionsRC$RC_END)
#if right side, then, (Reference sequence length - (Read Length - RC_END - insertion length)), 
#     if left side, it take the RC_END value. Track the parenthesis if the values are returning negative

            #What these lines just above are doing is correcting for previous scripts. They are taking into account that the    #
            #repeat motif that exists on the right of the break is actually the second instance of the repeated sequence.       #
            #This is subtracting the read length from the reference, giving a negative number, then adding back where previous  #
            #scripts thought the Direct Repeat Motifs stat/end are, plus the insertion length.                                  #

ConsistentInsertionsRC$motif_ID = paste(ConsistentInsertionsRC$ID, ConsistentInsertionsRC$RC_START_TRUE, ConsistentInsertionsRC$RC_END_TRUE, sep="-")
ConsistentInsertionsRC = ConsistentInsertionsRC[!duplicated(ConsistentInsertionsRC$motif_ID),]
            #The duplicates that get removed are due to SD-MMEJ events that could occur by either mechanism, therefore they     #
            #are included in both the direct repeat motifs and the reverse complement repeat motifs.                            #

ConsistentInsertionsRCR = subset(ConsistentInsertionsRC, SIDE == "RIGHT")
ConsistentInsertionsRCL = subset(ConsistentInsertionsRC, SIDE == "LEFT")

RCRPrimerPlotPoints = c() #Creates new data frome "RCRPrimerPlotPoints" for the for loop below
RCLPrimerPlotPoints = c() #Creates new data frome "RCLPrimerPlotPoints" for the for loop below

#----------------Reverse Complement - Left Side-----------------------------------------
if (nrow(ConsistentInsertionsRCL) > 0){
  for (i in 1:nrow(ConsistentInsertionsRCL)){
    ins <- as.character(ConsistentInsertionsRCL[i,35])    # create insertion sequence region to search
    #insertion
    dnains <- lapply(ins, DNAString)
    revins <- lapply(dnains, reverseComplement)
    ins <- lapply(revins, toString)
    rm <- toupper(ConsistentInsertionsRCL[i,43])
    #RM
    dr1 <- as.data.frame(str_locate_all(rm, as.character(ins))) # search deletion sequence for string
    if (nrow(dr1)!=1){
      dr2 <- dr1
      while (nrow(dr2)!=1){
        for (j in 1:nrow(dr1)){
          if (dr1[j,2]==nchar(as.character(ConsistentInsertionsRCL[i,43]))){ #RM
            break
          }else{
            prim <- paste(substring(ConsistentInsertionsRCL[i,43],1,dr1[j,1]-1), substring(ConsistentInsertionsRCL[i,43],dr1[j,2]+1,nchar(as.character(ConsistentInsertionsRCL[i,43]))),sep="")
            dnains2 <- lapply(prim, DNAString)
            revins2 <- lapply(dnains2, reverseComplement)
            prim <- lapply(revins2, toString)
            s <- dr1[j,1]-1
            der <- substring(as.character(ConsistentInsertionsRCL[i,34]), first=ConsistentInsertionsRCL[i,32]-(s-1), last=ConsistentInsertionsRCL[i,33]+(nchar(as.character(ConsistentInsertionsRCL[i,43]))-dr1[j,2]))
            #del_seq                              #left_del                          #right_del
            dr3 <- as.data.frame(str_locate_all(der, as.character(prim)))
            if (nrow(dr3)==0){ # if there is a match, then make dr3=dr2, which will make the nrow(dr2)=1, breaking the while loop
              next
            } else{
              dr2 <- dr3
            }
          }
        }
        break
      }
      RCLPrimerPlotPoints1 <- cbind(ConsistentInsertionsRCL[i,47], dr1[j,]) # once we know which one we want, we put that dr1 row
      #motif_ID
      RCLPrimerPlotPoints <- rbind(RCLPrimerPlotPoints,RCLPrimerPlotPoints1)
    } else {
      RCLPrimerPlotPoints1 <- cbind(ConsistentInsertionsRCL[i,47], dr1)
      #motif_ID
      RCLPrimerPlotPoints <- rbind(RCLPrimerPlotPoints,RCLPrimerPlotPoints1)
    }
  }
  
  
  colnames(RCLPrimerPlotPoints)[1] = "motif_ID"
  RCLPrimerPlotPoints$is_trans = "no"
  RCLPrimerPlotPoints$p1_start = (ConsistentInsertionsRCL$RC_END_TRUE)-(nchar(as.character(ConsistentInsertionsRCL$RM))-RCLPrimerPlotPoints$end-1)
  RCLPrimerPlotPoints$p1_end = ConsistentInsertionsRCL$RC_END_TRUE
  RCLPrimerPlotPoints$p1_length = 1+(RCLPrimerPlotPoints$p1_end - RCLPrimerPlotPoints$p1_start)
  RCLPrimerPlotPoints$mh_start = ConsistentInsertionsRCL$RC_START_TRUE
  RCLPrimerPlotPoints$mh_end = ConsistentInsertionsRCL$RC_START_TRUE+(RCLPrimerPlotPoints$start-2)
  RCLPrimerPlotPoints$mh_length = 1+(RCLPrimerPlotPoints$mh_end - RCLPrimerPlotPoints$mh_start)
  RCLPrimerPlotPoints$ins_start = ConsistentInsertionsRCL$RC_START_TRUE+(RCLPrimerPlotPoints$start-1)
  RCLPrimerPlotPoints$ins_end = ConsistentInsertionsRCL$RC_START_TRUE+(RCLPrimerPlotPoints$end-1)
  RCLPrimerPlotPoints$ins_length = 1+RCLPrimerPlotPoints$ins_end-RCLPrimerPlotPoints$ins_start
  RCLPrimerPlotPoints$p2_start = ConsistentInsertionsRCL$left_del - (RCLPrimerPlotPoints$p1_length-1)
  RCLPrimerPlotPoints$p2_end = ConsistentInsertionsRCL$left_del
  RCLPrimerPlotPoints$p2_length = 1+RCLPrimerPlotPoints$p2_end-RCLPrimerPlotPoints$p2_start
  ConsistentInsertionsRCL = merge(ConsistentInsertionsRCL,RCLPrimerPlotPoints, by="motif_ID")
  ConsistentInsertionsRCL$p_ID = paste(ConsistentInsertionsRCL$p1_start, ConsistentInsertionsRCL$p1_end,ConsistentInsertionsRCL$p2_start, ConsistentInsertionsRCL$p2_end, sep="-")
  RCLagg = aggregate(READS~p_ID, data=ConsistentInsertionsRCL, sum)
  RCLtable = as.data.frame(table(ConsistentInsertionsRCL$p_ID))
  RCLagg = merge(RCLagg, RCLtable, by.x="p_ID", by.y="Var1")
  ConsistentInsertionsRCL = merge(RCLagg, ConsistentInsertionsRCL, by="p_ID")
  ConsistentInsertionsRCL = subset(ConsistentInsertionsRCL, ! p1_start>p2_start & p1_start<p2_end)
  ConsistentInsertionsRCLsimplified = ConsistentInsertionsRCL[!duplicated(ConsistentInsertionsRCL$p_ID),]
  ConsistentInsertionsRCLsimplified$p2_start2 = ConsistentInsertionsRCLsimplified$p2_start-0.5
  ConsistentInsertionsRCLsimplified$p2_end2 = ConsistentInsertionsRCLsimplified$p2_end+0.5
  ConsistentInsertionsRCLsimplified = ConsistentInsertionsRCLsimplified[, -c(5:27), drop=FALSE]
  colnames(ConsistentInsertionsRCLsimplified)[2] = "READS"
  #Deletes some unnecessary columns to make the data frame easier to look through. Changes READS.x to READS   #
  #Theres more that can be removed still but leaving it at this point for now.                                #
}else{
  colnames=c('p_ID','READS','Freq','motif_ID','percent','percent_inaccurate','ID','DR_START','DR_END','RC_START','RC_END','consistency','left_del','right_del',
             'del_seq','insertion','plasmid','DRmotif_length','RCmotif_length','Loop.out','Snap.back','RIGHT_DEL','SIDE','RM','insertion_length','RC_START_TRUE','RC_END_TRUE',
             'start','end','is_trans','p1_start','p1_end','p1_length','mh_start','mh_end','mh_length','ins_start','ins_end',
             'ins_length','p2_start','p2_end','p2_length','p2_start2','p2_end2')
  ConsistentInsertionsRCLsimplified=data.frame(matrix(nrow=0, ncol=length(colnames)))
  colnames(ConsistentInsertionsRCLsimplified) = colnames
}


#----------------Reverse Complement Plot - Right Side-----------------------------------------
if (nrow(ConsistentInsertionsRCR) > 0){ 
  for (i in 1:nrow(ConsistentInsertionsRCR)){
    ins <- as.character(ConsistentInsertionsRCR[i,35])    # create insertion sequence region to search
    #insertion
    dnains <- lapply(ins, DNAString)
    revins <- lapply(dnains, reverseComplement)
    ins <- lapply(revins, toString)
    rm <- toupper(ConsistentInsertionsRCR[i,43]) #RM
    dr1 <- as.data.frame(str_locate_all(rm, as.character(ins))) # search deletion sequence for string
    if (nrow(dr1)!=1){
      dr2 <- dr1
      while (nrow(dr2)!=1){
        for (j in 1:nrow(dr1)){
          if (dr1[j,2]==nchar(as.character(ConsistentInsertionsRCR[i,43]))){ #RM
            break
          }else{
            prim <- paste(substring(ConsistentInsertionsRCR[i,43],1,dr1[j,1]-1), substring(ConsistentInsertionsRCR[i,43],dr1[j,2]+1,nchar(as.character(ConsistentInsertionsRCR[i,43]))),sep="")
            dnains2 <- lapply(prim, DNAString)
            revins2 <- lapply(dnains2, reverseComplement)
            prim <- lapply(revins2, toString)
            s <- dr1[j,1]-1
            der <- substring(as.character(ConsistentInsertionsRCR[i,34]), first=ConsistentInsertionsRCR[i,32]-(s-1), last=ConsistentInsertionsRCR[i,33]+(nchar(as.character(ConsistentInsertionsRCR[i,43]))-dr1[j,2]))
            #del_seq                              #left_del                               #right_del
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
      RCRPrimerPlotPoints1 <- cbind(ConsistentInsertionsRCR[i,47], dr1[k,]) # once we know which one we want, we put that dr1 row
      #motif_ID
      RCRPrimerPlotPoints <- rbind(RCRPrimerPlotPoints,RCRPrimerPlotPoints1)
    } else {
      RCRPrimerPlotPoints1 <- cbind(ConsistentInsertionsRCR[i,47], dr1)
      #motif_ID
      RCRPrimerPlotPoints <- rbind(RCRPrimerPlotPoints,RCRPrimerPlotPoints1)
    }
  }
  
  colnames(RCRPrimerPlotPoints)[1] = "motif_ID"
  RCRPrimerPlotPoints$is_trans = "no"
  RCRPrimerPlotPoints$p1_start = ConsistentInsertionsRCR$RC_START_TRUE
  RCRPrimerPlotPoints$p1_end = ConsistentInsertionsRCR$RC_START_TRUE+(RCRPrimerPlotPoints$start-2)
  RCRPrimerPlotPoints$p1_length = 1+RCRPrimerPlotPoints$p1_end-RCRPrimerPlotPoints$p1_start
  RCRPrimerPlotPoints$mh_start = ConsistentInsertionsRCR$RC_END_TRUE-(nchar(as.character(ConsistentInsertionsRCR$RM))-RCRPrimerPlotPoints$end-1)
  RCRPrimerPlotPoints$mh_end = ConsistentInsertionsRCR$RC_END_TRUE
  RCRPrimerPlotPoints$mh_length = 1+ RCRPrimerPlotPoints$mh_end-RCRPrimerPlotPoints$mh_start
  RCRPrimerPlotPoints$ins_start = ConsistentInsertionsRCR$RC_START_TRUE+RCRPrimerPlotPoints$start-1
  RCRPrimerPlotPoints$ins_end = ConsistentInsertionsRCR$RC_START_TRUE+(RCRPrimerPlotPoints$end-1)
  RCRPrimerPlotPoints$ins_length = 1+RCRPrimerPlotPoints$ins_end-RCRPrimerPlotPoints$ins_start
  RCRPrimerPlotPoints$p2_start = ConsistentInsertionsRCR$RIGHT_DEL+1
  RCRPrimerPlotPoints$p2_end = (ConsistentInsertionsRCR$RIGHT_DEL+1)+(RCRPrimerPlotPoints$p1_length-1)
  RCRPrimerPlotPoints$p2_length = 1+RCRPrimerPlotPoints$p2_end-RCRPrimerPlotPoints$p2_start
  ConsistentInsertionsRCR = merge(ConsistentInsertionsRCR, RCRPrimerPlotPoints, by="motif_ID")
  ConsistentInsertionsRCR$p_ID = paste(ConsistentInsertionsRCR$p1_start, ConsistentInsertionsRCR$p1_end ,ConsistentInsertionsRCR$p2_start, ConsistentInsertionsRCR$p2_end, sep="-")
  RCRagg = aggregate(READS~p_ID, data=ConsistentInsertionsRCR, sum)
  RCRtable = as.data.frame(table(ConsistentInsertionsRCR$p_ID))
  RCRagg = merge(RCRagg, RCRtable, by.x="p_ID", by.y="Var1")
  ConsistentInsertionsRCR = merge(RCRagg, ConsistentInsertionsRCR, by="p_ID")
  ConsistentInsertionsRCR = subset(ConsistentInsertionsRCR, ! p2_start>p1_start & p2_start<p1_end)
  ConsistentInsertionsRCRsimplified = ConsistentInsertionsRCR[!duplicated(ConsistentInsertionsRCR$p_ID),]
  ConsistentInsertionsRCRsimplified$p2_start2 = ConsistentInsertionsRCRsimplified$p2_start-0.5
  ConsistentInsertionsRCRsimplified$p2_end2 = ConsistentInsertionsRCRsimplified$p2_end+0.5
  ConsistentInsertionsRCRsimplified = ConsistentInsertionsRCRsimplified[, -c(5:27), drop=FALSE]
  colnames(ConsistentInsertionsRCRsimplified)[2] = "READS"
  
  #Deletes some unnecessary columns to make the data frame easier to look through. Changes READS.x to READS   #
              #Theres more that can be removed still but leaving it at this point for now.                                #
}else{
  colnames=colnames(ConsistentInsertionsRCLsimplified)
  ConsistentInsertionsRCRsimplified = data.frame(matrix(nrow=0, ncol=length(colnames)))
  colnames(ConsistentInsertionsRCRsimplified) = colnames
}


#----------------Combining Direct Repeat and Reverse Complement Repeat-----------------------------------------

DirectRepeats = rbind(ConsistentInsertionsDRLsimplified, ConsistentInsertionsDRRsimplified)

try( {
  DirectRepeats$RC_START_TRUE = "NA"
  DirectRepeats$RC_END_TRUE = "NA"
  DirectRepeats$mechanism = ifelse(DirectRepeats$is_trans=="trans", "Trans", "Loop-out")
}, silent=TRUE)


ReverseCompRepeats = rbind(ConsistentInsertionsRCLsimplified, ConsistentInsertionsRCRsimplified)
try( {
  ReverseCompRepeats$DR_START_TRUE = "NA"
  ReverseCompRepeats$DR_END_TRUE = "NA"
  ReverseCompRepeats$mechanism = "Snap-back"
}, silent=TRUE)


AllRepeats = rbind(ReverseCompRepeats, DirectRepeats)

AllRepeats$p_mechID = paste(AllRepeats$p_ID,AllRepeats$mechanism,sep="-")
        # At this point since we already aggregated the reads by IDs and simplified by removing duplicate IDs, each entry     #
        # should have its own p2_mechID. If there are duplicates for some reason, run the 5 lines below that have # in front. #

#AllRepeatsAggr = aggregate(READS~p2_mechID,data=AllRepeats, sum)
#AllRepeatsTable=as.data.frame(table(AllRepeats$p2_mechID))
#AllRepeatsAggrTable = merge(AllRepeatsAggr, AllRepeatsTable, by.x="p2_mechID", by.y = "Var1")
#AllRepeats = merge(AllRepeatsAggrTable, AllRepeats, by="p2_mechID")
#AllRepeatsSimplified = AllRepeats[!duplicated(AllRepeats$p2_mechID),]

AllRepeats$p1_start2 <- AllRepeats$p1_start-.5
AllRepeats$p1_end2 <- AllRepeats$p1_end+.5
AllRepeats$mh_start2 <- AllRepeats$mh_start-.5
AllRepeats$mh_end2 <- AllRepeats$mh_end+.5
AllRepeats$ins_start2 <- AllRepeats$ins_start-.5
AllRepeats$ins_end2 <- AllRepeats$ins_end+.5
AllRepeats$percent_insertion_jxn <- AllRepeats$READS/sum(AllRepeats$READS)*100

#----------------Single-step Insertion Reads Primer Plot Data Manipulation-------------------

if (nrow(AllRepeats)>0){
  AllRepeats$jxnID = paste(AllRepeats$p1_start, AllRepeats$p1_end, AllRepeats$p2_start, AllRepeats$p2_end, AllRepeats$mechanism, sep="-")
  # At this point since we already aggregated the reads by IDs and simplified by removing duplicate IDs, each entry   #
  # should have its own jxnID. If there are duplicates for some reason, run the 5 lines below that have # in front.   #
  
  #jxnAGGR =  aggregate(READS~jxnID, data=AllRepeats, sum)
  #jxnTable = as.data.frame(table(AllRepeats$jxnID))
  #jxnAGGRTable = merge(jxnAGGR, jxnTable, by.x="jxnID", by.y="Var1")
  #AllRepeats = merge(jxnAGGRTable, AllRepeats, by="jxnID")
  #jxnSimplified = AllRepeats[!duplicated(AllRepeats$jxnID),]
  
  jxnLeft = subset(AllRepeats, SIDE=="LEFT")
  jxnRight = subset(AllRepeats, SIDE=="RIGHT")
  jxnLeft = jxnLeft[order(jxnLeft$p2_start),]
  jxnRight = jxnRight[order(jxnRight$p2_start, decreasing=TRUE),]
  
  # can be empty
  try({jxnLeft$order = paste(10:(9+nrow(jxnLeft)), "-", sep = "")}, silent=T)
  try({jxnRight$order = paste(10:(9+nrow(jxnRight)), "-", sep = "")}, silent=T)
  
  jxn = rbind(jxnLeft, jxnRight)
  try({jxn = jxn[order(jxn$order, decreasing=TRUE),]}, silent=T)
  jxn = arrange(transform(jxn, jxnID=factor(jxnID, levels = jxnID)),jxnID)
  
  ## Remove not actual insertion events, most of these are complex events that are actually simple deletions.
  ## Nick will send Rebecca some examples to check where they are being misidentified
  jxn3 = subset(jxn, p1_start2!=p2_start2)
  jxn3=subset(jxn3, abs(p1_start-p2_start)>=2)
  jxn3=subset(jxn3, abs(p1_start-p2_end)>2)
  jxn3=subset(jxn3, p1_start2!=p1_end2)
  jxn3$percent_insertion_jxn <- jxn3$READS/sum(jxn3$READS)*100
  write.csv (jxn3, paste(outdir, "/", "table_outputs/", plasmid, "_break.csv", sep=""))
  
  #----------------Single-step Insertion Reads Primer Plot-------------------
  rect_left <- c(33.5,43.5,53.5,63.5,73.5,83.5,93.5,103.5,113.5,123.5,133.5,143.5,153.5,163.5,173.5,183.5,193.5)
  rectangles <- data.frame(xmin = rect_left, xmax = rect_left + 5, ymin = -Inf,ymax = Inf)
  PrimerPlotBreaks = min(jxn3$p1_start2-2.5): max(jxn3$p1_end2+2.5)
  PrimerPlotLabel = dput(as.character(referenceSplit[PrimerPlotBreaks]))
  PrimerPlot <- ggplot()+geom_boxplot(data=jxn3, aes(x=jxnID, ymin = p2_start2, lower = p2_start2, middle = p2_start2, 
                                                     upper = p2_end2, ymax = p1_start,fill=percent_insertion_jxn, linetype=mechanism),
                                      colour="black",stat = "identity",width=.8) + coord_flip()+scale_linetype_manual(values=c("solid", "dashed","dotted"),guide="none")+
    scale_fill_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),values=c(1.0,0.8,0.6,0.4,0.2,0),
                         limits=c(0,max(jxn3$percent_insertion_jxn+1.5)),guide_legend(title="% Single-step\nInsertion Reads"))+
    scale_y_continuous(limits=c(min(jxn3$p1_start2-2.5),max(jxn3$p1_end2+2.5)),breaks=PrimerPlotBreaks,
                       labels = PrimerPlotLabel,
                       name = "Nucleotide",expand=c(0,0))+theme_classic(base_size = 10)+theme(axis.text.y= element_blank(),axis.title.y=element_blank(),
                                                                                              axis.ticks=element_blank())+geom_rect(data=rectangles, aes(ymin=as.numeric(xmin), ymax=as.numeric(xmax), xmin=ymin, xmax=ymax),
                                                                                                                                    fill="grey", alpha=0.4, colour=NA)+geom_boxplot(data=jxn3, aes(x=jxnID, ymin = p2_start2, lower = p2_start2, middle = p2_start2, 
                                                                                                                                                                                                   upper = p2_end2, ymax = p2_end2,fill=percent_insertion_jxn, linetype=mechanism),colour="black",stat = "identity",width=.8)+ 
    geom_boxplot(data=jxn3, aes(x=jxnID, ymin = p1_start2, lower = p1_start2, middle = p1_start2, upper = p1_end2, ymax = p1_end2,
                                fill=percent_insertion_jxn, linetype=mechanism),colour="black",stat = "identity",width=.8)+geom_hline(yintercept = BreakPointFromLeft+0.5, colour="red", size=0.6)
  PrimerPlot
  pdf(paste(outdir, "/", "plots/", plasmid, "_Ins_Primer_Plot.pdf", sep=""), height = 5, width = 10)
  PrimerPlot
  dev.off()

  #----------------Insertion Repeat Motif Plot - Data Manipulation-------------------
  DRMotif = rbind(ConsistentInsertionsDRL, ConsistentInsertionsDRR)
  RCMotif = rbind(ConsistentInsertionsRCL, ConsistentInsertionsRCR)
  
  if (nrow(DRMotif)>0){
    DRMotif$RC_START_TRUE = "NA"
    DRMotif$RC_END_TRUE = "NA"
    DRMotif$mechanism = ifelse(DRMotif$is_trans=="trans", "Trans", "Loop-out")
  }
  if (nrow(RCMotif)>0){
    RCMotif$DR_START_TRUE = "NA"
    RCMotif$DR_END_TRUE = "NA"
    RCMotif$mechanism = "Snap-back"
  }
  InsRepeatMotif = rbind(DRMotif, RCMotif)
  InsRepeatMotif$MOTIF_START <- ifelse(InsRepeatMotif$DR_START_TRUE=="NA", InsRepeatMotif$RC_START_TRUE,
                                       InsRepeatMotif$DR_START_TRUE)
  InsRepeatMotif$MOTIF_END <- ifelse(InsRepeatMotif$DR_END_TRUE=="NA",InsRepeatMotif$RC_END_TRUE,
                                     InsRepeatMotif$DR_END_TRUE)
  InsRepeatMotif$motif <- paste(InsRepeatMotif$MOTIF_START,InsRepeatMotif$MOTIF_END,sep="-")
  InsRepeatMotif$motif_mechID <- paste(InsRepeatMotif$motif,InsRepeatMotif$mechanism,sep="-")
  write.csv(InsRepeatMotif,paste(outdir, "/", "table_outputs/", plasmid, "_insertion_repeat_motif.csv", sep=""))
  
  InsRepeatMotifLeft = subset(InsRepeatMotif, SIDE=="LEFT")
  InsRepeatMotifRight = subset(InsRepeatMotif, SIDE=="RIGHT")
  InsRepeatMotifLeft = InsRepeatMotifLeft[order(as.numeric(InsRepeatMotifLeft$MOTIF_START), decreasing = TRUE),]
  InsRepeatMotifRight = InsRepeatMotifRight[order(as.numeric(InsRepeatMotifRight$MOTIF_START), decreasing = TRUE),]
  
  InsRepeatMotifLeft$order = paste(10:(9+nrow(InsRepeatMotifLeft)), "-", sep="")
  InsRepeatMotifRight$order = paste(10:(9+nrow(InsRepeatMotifRight)), "-", sep="")
  
  InsRepeatMotif = rbind(InsRepeatMotifLeft, InsRepeatMotifRight)
  InsRepeatMotif = InsRepeatMotif[order(InsRepeatMotif$order, decreasing = TRUE),]
  InsRepeatMotif$MOTIF_START2 = as.numeric(InsRepeatMotif$MOTIF_START)-0.5
  InsRepeatMotif$MOTIF_END2 = as.numeric(InsRepeatMotif$MOTIF_END)+0.5
  
  ## multiple repeat motifs can have same primer pairs
  ## this code is also removing complex that are actually simple deletions
  InsRepeatMotif = subset(InsRepeatMotif, p1_start!=p2_start)
  InsRepeatMotif = subset(InsRepeatMotif, !(p2_end>p1_start & p2_end<p1_end))
  InsRepeatMotif = subset(InsRepeatMotif, abs(p1_start-p2_start)>=2)
  InsRepeatMotif = subset(InsRepeatMotif, abs(p1_start-p2_end)>2)
  InsRepeatMotif1 = subset(InsRepeatMotif, MOTIF_END2>BreakPointFromLeft+0.5 & MOTIF_START2<BreakPointFromLeft+0.5)
  InsRepeatMotif = setdiff(InsRepeatMotif, InsRepeatMotif1)
  InsRepeatMotif$percent_insertion = (InsRepeatMotif$READS.y/sum(InsRepeatMotif$READS.y))*100
  
  #----------------Insertion Repeat Motif Plot-------------------
  rect_left <- c(43.5,53.5,63.5,73.5,83.5,93.5,103.5,113.5,123.5,133.5,143.5,153.5,163.5,173.5,183.5,193.5,203.5)
  rectangles <- data.frame(xmin = rect_left,xmax = rect_left + 5,ymin = -Inf,ymax = Inf)
  InsRepeatMotifPlotBreaks = min(InsRepeatMotif$MOTIF_START2-1.5):max(InsRepeatMotif$MOTIF_END2+1.5)
  InsRepeatMotifPlotLabel = dput(as.character(referenceSplit[InsRepeatMotifPlotBreaks]))
  RepeatMotifPlot <- ggplot()+
    geom_boxplot(data=InsRepeatMotif, aes(x=motif_ID, ymin = MOTIF_START2, lower = MOTIF_START2, middle = MOTIF_START2, 
                                          upper = MOTIF_END2, ymax = MOTIF_END2,fill=percent_insertion, linetype=mechanism),
                 colour="white", stat = "identity",width=.8,lwd=.5) + coord_flip()+scale_linetype_manual(values=c("solid", "dashed","dotted"),guide="none")+
    scale_fill_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),values=c(1.0,0.8,0.6,0.4,0.2,0),
                         guide_legend(title="% Single-step\nInsertion Reads"),limits=c(0,max(InsRepeatMotif$percent_insertion+0.5)))+scale_y_continuous(limits=c(min(InsRepeatMotif$MOTIF_START2-1.5),max(InsRepeatMotif$MOTIF_END2+1.5)),
                                                                                                                                                        breaks=InsRepeatMotifPlotBreaks,
                                                                                                                                                        labels = InsRepeatMotifPlotLabel,
                                                                                                                                                        name = "Nucleotide",expand=c(0,0))+theme_classic(base_size = 10)+theme(legend.text=element_text(size=8,face="bold"),
                                                                                                                                                                                                                               axis.text.y= element_blank(),axis.title.y=element_blank(),axis.text.x= element_text(face="bold"),axis.title.x=element_text(face="bold"),
                                                                                                                                                                                                                               axis.ticks=element_blank(),
                                                                                                                                                                                                                               panel.border = element_rect(colour = "black", fill=NA, size=1))+guides(guide_legend(title.theme = element_text(size=8, face="bold", angle=0),
                                                                                                                                                                                                                                                                                                                   label.theme=element_text(size=8, face="bold", angle=0)))+geom_rect(data=rectangles, aes(ymin=as.numeric(xmin), 
                                                                                                                                                                                                                                                                                                                                                                                                           ymax=as.numeric(xmax), xmin=ymin, xmax=ymax),fill="grey", alpha=0.4, colour=NA)+
    geom_boxplot(data=InsRepeatMotif, aes(x=motif_ID, ymin = MOTIF_START2, lower = MOTIF_START2, middle = MOTIF_START2, 
                                          upper = MOTIF_END2, ymax = MOTIF_END2,fill=percent_insertion, linetype=mechanism), colour="black",stat = "identity",
                 width=.8,lwd=.5)+ geom_hline(yintercept = BreakPointFromLeft+0.5, colour="red", size=0.75)
  RepeatMotifPlot
  pdf(paste(outdir, "/", "plots/", plasmid, "_Insertion_Repeat_Motif_Plot.pdf", sep=""), height = 5, width = 10)
  RepeatMotifPlot
  dev.off()
  
  #----------------Insertion Resection - Data Manipulation-----------
  #To generate a plot of the furthest identifiable point of resection or duplex unwinding = furthest from break repeat motif
  InsLeft = subset(InsRepeatMotif, SIDE=="LEFT")
  InsLeftLoop= subset(InsLeft, mechanism=="Loop-out")
  InsLeftSnap = subset(InsLeft, mechanism=="Snap-back")
  InsLeftLoop$RM_to_break = BreakPointFromLeft-as.numeric(InsLeftLoop$p1_start)
  InsLeftSnap$RM_to_break = BreakPointFromLeft-as.numeric(InsLeftSnap$p1_start)
  InsLeft=rbind(InsLeftLoop, InsLeftSnap)
  InsRight = subset(InsRepeatMotif, SIDE=="RIGHT")
  InsRightLoop = subset(InsRight, mechanism=="Loop-out")
  InsRightSnap = subset(InsRight, mechanism=="Snap-back")
  InsRightLoop$RM_to_break = as.numeric(InsRightLoop$p1_end)-BreakPointFromLeft
  InsRightSnap$RM_to_break = as.numeric(InsRightSnap$p1_end)-BreakPointFromLeft
  InsRight = rbind(InsRightLoop, InsRightSnap)
  InsFromBreak = rbind(InsLeft, InsRight)
  InsFromBreak = InsFromBreak[,-c(2:4, 6:19, 23:26, 28, 29, 51:65)]
  InsFromBreak$percent= (InsFromBreak$READS.y/sum(InsFromBreak$READS.y))*100
  write.csv(InsFromBreak, paste(outdir, "/", "table_outputs/", plasmid, "_insertion_resection_data.csv", sep=""))
  InsFromBreakMean=sum(InsFromBreak$RM_to_break*InsFromBreak$READS.y)/sum(InsFromBreak$READS.y)
  
  #----------------Insertion Resection - Plot - Side-----------
  InsertionResectionPlotSide = ggplot(InsFromBreak)+
    geom_bar(aes(x=RM_to_break, y=percent, color=SIDE), position="stack", stat="identity", fill="grey50") +
    theme_bw()+
    geom_vline(aes(xintercept = InsFromBreakMean), colour="blue", size=0.75,linetype = "longdash")+
    #geom_vline(aes(xintercept = median), colour="red", size=0.75,linetype = "longdash",data=distMedian)+
    scale_y_continuous(name = "Percent Insertions (%)")+  
    scale_x_continuous(name="Resection/Unwinding Required for SD-MMEJ Action (bp)",
                       breaks = seq(0, max(InsFromBreak$RM_to_break), by = 2))+
    theme(axis.text.x=element_text(size=10, face="bold"),
          axis.text.y=element_text(size=10, face="bold"),
          axis.title.x=element_text(size=12, face="bold"),
          axis.title.y=element_text(size=12, face="bold"))+
    theme(strip.text.x = element_text(size=10, face="bold"),
          strip.text.y = element_text(size=10, face="bold"))
  InsertionResectionPlotSide
  pdf(paste(outdir, "/", "plots/", plasmid, "_Insertion_Resection_Plot_Side.pdf", sep=""), height = 5, width = 10)
  InsertionResectionPlotSide
  dev.off()
  
  #----------------Insertion Resections - Plot - Mechanism-----------
  InsertionResectionPlotMechanism = ggplot(InsFromBreak)+
    geom_bar(aes(x=RM_to_break, y=percent, color=mechanism), position="stack", stat="identity", fill="grey50") +
    theme_bw()+
    geom_vline(aes(xintercept = InsFromBreakMean), colour="blue", size=0.75,linetype = "longdash")+
    #geom_vline(aes(xintercept = median), colour="red", size=0.75,linetype = "longdash",data=distMedian)+
    scale_y_continuous(name = "Percent Insertions (%)")+  
    scale_x_continuous(name="Resection/Unwinding Required for SD-MMEJ Action (bp)",
                       breaks = seq(0, max(InsFromBreak$RM_to_break), by = 2))+
    theme(axis.text.x=element_text(size=10, face="bold"),
          axis.text.y=element_text(size=10, face="bold"),
          axis.title.x=element_text(size=12, face="bold"),
          axis.title.y=element_text(size=12, face="bold"))+
    theme(strip.text.x = element_text(size=10, face="bold"),
          strip.text.y = element_text(size=10, face="bold"))
  InsertionResectionPlotMechanism
  pdf(paste(outdir, "/", "plots/", plasmid, "_Insertion_Resection_Plot_Mechanism.pdf", sep=""), height = 5, width = 10)
  InsertionResectionPlotMechanism
  dev.off()
  
  #----------------Insertion Flap Plot - Data Manipulation-----------------------------------
  InsFlap=InsRepeatMotif
  InsFlap$LeftFlap = InsFlap$left_del+0.5
  InsFlap$RightFlap = InsFlap$RIGHT_DEL+0.5
  InsFlap$FlapID = paste(InsFlap$LeftFlap, InsFlap$RightFlap, InsFlap$mechanism)
  InsFlapAggr=aggregate(READS.y~FlapID, data=InsFlap,sum)
  InsFlapTable=as.data.frame(table(InsFlap$FlapID))
  InsFlapAggr=merge(InsFlapAggr, InsFlapTable, by.x="FlapID", by.y="Var1")
  InsFlap=merge(InsFlapAggr, InsFlap, by="FlapID")
  InsFlap=InsFlap[!duplicated(InsFlap$FlapID),]
  colnames(InsFlap)[2]= "READS"
  InsFlap$percent_flap = InsFlap$READS/(sum(InsFlap$READS))*100
  
  #----------------Insetion Flap Plot - Plot-----------------------------------
  rect_left <- c(43.5,53.5,63.5,73.5,83.5,93.5,103.5,113.5,123.5,133.5,143.5,153.5,163.5,173.5,183.5,193.5,203.5)
  InsFlapPlotBreaks = min(InsFlap$LeftFlap-1.5):max(InsFlap$RightFlap+1.5)
  InsFlapPlotLabel = dput(as.character(referenceSplit[InsFlapPlotBreaks]))
  InsFlapPlot <- ggplot()+
    geom_boxplot(data=InsFlap, aes(x=motif_ID, ymin = LeftFlap, lower = LeftFlap, middle = LeftFlap, 
                                   upper = RightFlap, ymax = RightFlap,fill=percent_flap, linetype=mechanism),
                 colour="white", stat = "identity",width=.8,lwd=.5) + coord_flip()+scale_linetype_manual(values=c("solid", "dashed","dotted"),guide="none")+
    scale_fill_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),values=c(1.0,0.8,0.6,0.4,0.2,0),
                         guide_legend(title="% Deletion From\nSingle Step\nInsertion Events"),limits=c(0,max(InsFlap$percent_flap+1.5)))+scale_y_continuous(limits=c(min(InsFlap$LeftFlap-1.5),max(InsFlap$RightFlap+1.5)),
                                                                                                                                                            breaks=InsFlapPlotBreaks,
                                                                                                                                                            labels = InsFlapPlotLabel,
                                                                                                                                                            name = "Nucleotide",expand=c(0,0))+theme_classic(base_size = 10)+theme(legend.text=element_text(size=8,face="bold"),
                                                                                                                                                                                                                                   axis.text.y= element_blank(),axis.title.y=element_blank(),axis.text.x= element_text(face="bold"),axis.title.x=element_text(face="bold"),
                                                                                                                                                                                                                                   axis.ticks=element_blank(),
                                                                                                                                                                                                                                   panel.border = element_rect(colour = "black", fill=NA, size=1))+guides(guide_legend(title.theme = element_text(size=8, face="bold", angle=0),
                                                                                                                                                                                                                                                                                                                       label.theme=element_text(size=8, face="bold", angle=0)))+geom_rect(data=rectangles, aes(ymin=as.numeric(xmin), 
                                                                                                                                                                                                                                                                                                                                                                                                               ymax=as.numeric(xmax), xmin=ymin, xmax=ymax),fill="grey", alpha=0.4, colour=NA)+
    geom_boxplot(data=InsFlap, aes(x=motif_ID, ymin = LeftFlap, lower = LeftFlap, middle = LeftFlap, 
                                   upper = RightFlap, ymax = RightFlap,fill=percent_flap, linetype=mechanism), colour="black",stat = "identity",
                 width=.8,lwd=.5)+ geom_hline(yintercept = BreakPointFromLeft+0.5, colour="red", size=0.75)
  
  InsFlapPlot
  pdf(paste(outdir, "/", "plots/", plasmid, "_Insertion_Flap_Plot.pdf", sep=""), height = 5, width = 10)
  InsFlapPlot
  dev.off()
  
  #----------------Insertion Repeat Motif by Side - Plot---------------
  InsSide = InsRepeatMotif
  InsSidePlot = ggplot(InsSide)+
    geom_bar(aes(x=SIDE, y=percent_insertion, color=mechanism), position="stack", stat="identity", fill="grey50") +
    theme_bw()+
    scale_y_continuous(name = "Percent Usage (%)")+  
    scale_x_discrete(name="Side")+
    theme(axis.text.x=element_text(size=10, face="bold"),
          axis.text.y=element_text(size=10, face="bold"),
          axis.title.x=element_text(size=12, face="bold"),
          axis.title.y=element_text(size=12, face="bold"))+
    theme(strip.text.x = element_text(size=10, face="bold"),
          strip.text.y = element_text(size=10, face="bold"))
  InsSidePlot
  pdf(paste(outdir, "/", "plots/", plasmid, "_Insertion_Side_Usage_Plot.pdf", sep=""), height = 5, width = 10)
  InsSidePlot
  dev.off()
  
  ## Second set of insertion plots
  #--------------------Primer Distance and Length - Data Manipulation--------------------------
  # csv input - _break.csv exported at end of "Combining Direct Repeat and Reverse Complement Repeat" Section
  #InsBreaks = read.csv(paste(outdir, "/", "table_outputs/", plasmid, "_break.csv", sep=""), row.names = 1)
  InsBreaks = jxn3
  InsBreaks$PG = plasmid
  no_trans = subset(InsBreaks, is_trans=="no")
  #trans = subset(InsBreaks, is_trans=="trans")
  right = subset(no_trans, SIDE=="RIGHT")
  left = subset(no_trans, SIDE=="LEFT")
  right$p1_p2 = right$p1_start-(right$p2_end+1)
  left$p1_p2 = left$p2_start-(left$p1_end+1)
  #trans$p1_p2 = 0
  #dist = rbind(right,left,trans)
  dist = rbind(right,left)
  dist.nt = subset(dist, mechanism!="Trans")
  dist.nt = subset(dist.nt, !(is.na(dist.nt["p1_length"])))
  lengthMean = data.frame(plasmid=plasmid, mean(dist.nt$p1_length))
  colnames(lengthMean)[2]="mean"
  lengthMedian = data.frame(plasmid=plasmid, median(dist.nt$p1_length))
  colnames(lengthMedian)[2]="median"
  dist.nt$p1_p2_abs = abs(dist.nt$p1_p2)
  distMean = data.frame(plasmid=plasmid, mean(dist.nt$p1_p2_abs))
  colnames(distMean)[2] = "mean"
  distMedian = data.frame(plasmid=plasmid, median(dist.nt$p1_p2_abs))
  colnames(distMedian)[2] = "median"
  len.mean = data.frame(plasmid=plasmid, mean(dist.nt$mh_length))
  colnames(len.mean)[2] = "mean"
  len.median = data.frame(plasmid=plasmid, median(dist.nt$mh_length))
  colnames(len.median)[2] = "median"
  
  PrimerData = jxn3[,c(1,2,33,31,32,40,41,23)]
  colnames(PrimerData)[3]="Primer_Length"
  PrimerDataLeft=subset(PrimerData, SIDE=="LEFT")
  PrimerDataRight=subset(PrimerData, SIDE=="RIGHT")
  PrimerDataLeft$Distance=PrimerDataLeft$p2_start-1-PrimerDataLeft$p1_end
  PrimerDataRight$Distance=PrimerDataRight$p1_start-1-PrimerDataRight$p2_end
  PrimerData=rbind(PrimerDataLeft, PrimerDataRight)
  
  PrimerLength = aggregate(READS~Primer_Length, data=PrimerData, sum)
  PrimerLength$percent=PrimerLength$READS/sum(PrimerLength$READS)*100
  write.csv(PrimerLength, paste(outdir, "/", "table_outputs/", plasmid, "_Primer_Length.csv", sep=""))
  
  PrimerDistance = aggregate(READS~Distance, data=PrimerData, sum)
  PrimerDistance$percent = PrimerDistance$READS/sum(PrimerDistance$READS)*100
  write.csv(PrimerDistance, paste(outdir, "/", "table_outputs/", plasmid, "_Primer_Distance.csv", sep=""))
  
  #--------------------Primer Distance Plot - Inaccurate Reads--------------------------
  PrimerDistancePlotInaccurateReads = ggplot(dist.nt)+
    geom_bar(aes(x=p1_p2_abs, y=as.numeric(percent_inaccurate), color=mechanism), position="stack", stat="identity", fill="grey50") +
    theme_bw()+
    #facet_wrap(~PG, scales="free_y")+
    geom_vline(aes(xintercept = mean), colour="blue", size=0.75,linetype = "longdash",data=distMean)+
    geom_vline(aes(xintercept = median), colour="red", size=0.75,linetype = "longdash",data=distMedian)+
    scale_y_continuous(name = "Percent Inaccurate Reads")+  
    scale_x_continuous(name="Distance between Primer Pairs (bp)",
                       breaks = seq(0, 60, by = 2))+
    theme(axis.text.x=element_text(size=10, face="bold"),
          axis.text.y=element_text(size=10, face="bold"),
          axis.title.x=element_text(size=12, face="bold"),
          axis.title.y=element_text(size=12, face="bold"))+
    theme(strip.text.x = element_text(size=10, face="bold"),
          strip.text.y = element_text(size=10, face="bold"))
  PrimerDistancePlotInaccurateReads
  pdf(paste(outdir, "/", "plots/", plasmid, "_Primer_Distance_Plot_Inaccurate_Reads.pdf", sep=""), width=15)
  PrimerDistancePlotInaccurateReads
  dev.off()
  
  #--------------------Primer Distance Plot - Insertion Events--------------------------
  PrimerDistancePlotInsertion = ggplot(dist.nt)+
    geom_bar(aes(x=p1_p2_abs, y=percent_insertion_jxn, color=mechanism), position="stack", stat="identity", fill="grey50") +
    theme_bw()+
    #facet_wrap(~PG, scales="free_y")+
    geom_vline(aes(xintercept = mean), colour="blue", size=0.75,linetype = "longdash",data=distMean)+
    geom_vline(aes(xintercept = median), colour="red", size=0.75,linetype = "longdash",data=distMedian)+
    scale_y_continuous(name = "Percent Insertion Events")+  
    scale_x_continuous(name="Distance between Primer Pairs (bp)",
                       breaks = seq(0, 60, by = 2))+
    theme(axis.text.x=element_text(size=10, face="bold"),
          axis.text.y=element_text(size=10, face="bold"),
          axis.title.x=element_text(size=12, face="bold"),
          axis.title.y=element_text(size=12, face="bold"))+
    theme(strip.text.x = element_text(size=10, face="bold"),
          strip.text.y = element_text(size=10, face="bold"))
  PrimerDistancePlotInsertion
  pdf(paste(outdir, "/", "plots/", plasmid, "_Primer_Distance_Plot_Insertion.pdf", sep=""), width=15)
  PrimerDistancePlotInsertion
  dev.off()
  
  #--------------------Primer Length Plots - Insertion Events--------------------------
  PrimerLengthPlotInsertion <- ggplot(dist.nt)+
    geom_bar(aes(x=p2_length, y=percent_insertion_jxn, color=mechanism), position="stack", stat="identity", fill="grey50") +
    theme_bw()+
    #facet_wrap(~PG, scales="free_y")+
    geom_vline(aes(xintercept = mean), colour="blue", size=0.75,linetype = "longdash",data=lengthMean)+
    geom_vline(aes(xintercept = median), colour="red", size=0.75,linetype = "longdash",data=lengthMedian)+
    scale_y_continuous(name = "Percent Insertion Events")+  
    scale_x_continuous(name="Primer Length (bp)",
                       breaks = seq(1, 20, by = 1))+
    theme(axis.text.x=element_text(size=10, face="bold"),
          axis.text.y=element_text(size=10, face="bold"),
          axis.title.x=element_text(size=12, face="bold"),
          axis.title.y=element_text(size=12, face="bold"))+
    theme(strip.text.x = element_text(size=10, face="bold"),
          strip.text.y = element_text(size=10, face="bold"))
  PrimerLengthPlotInsertion
  pdf(paste(outdir, "/", "plots/", plasmid, "_Primer_Length_Plot_Insertion.pdf", sep=""), width=15)
  PrimerLengthPlotInsertion
  dev.off()
  
  ## Third set of insertion plots
  
  #--------------------------Insertion Length Plot - All----------------
  AllInsertions = subset(all, CLASS == "Insertion")
  ReadsAggInsertionLength = aggregate(READS~INSERTION_LENGTH + CONSISTENCY, data=AllInsertions, sum)
  ReadsAggInsertionLength$Percent_Insertion = ReadsAggInsertionLength$READS/ sum(ReadsAggInsertionLength$READS)*100
  ReadsAggInsertionLength$READSxLength = ReadsAggInsertionLength$READS * ReadsAggInsertionLength$INSERTION_LENGTH
  AllInsLenMean = sum(ReadsAggInsertionLength$READSxLength)/sum(ReadsAggInsertionLength$READS)
  
  ReadsAggInsertionLength = ReadsAggInsertionLength[-c(24),]
  
  len.median = data.frame(plasmid="Iw7_Flex", median(dist.nt$mh_length))
  
  AllInsertionLengthPlot = ggplot(ReadsAggInsertionLength) + geom_bar(aes(x=INSERTION_LENGTH, y=Percent_Insertion, fill= CONSISTENCY), stat="identity", colour="black")+
    theme_bw()+scale_y_continuous(name = "Percent of Insertions")+ scale_x_continuous(name = "Insertion Length (bp)", breaks = seq(0, max(ReadsAggInsertionLength$INSERTION_LENGTH+1), by = 2))+
    #facet_wrap(~"WT", scales="free_y")+
    theme(plot.title = element_text(color="grey28"),
          axis.text.x=element_text(size=10, face="bold"),
          axis.text.y=element_text(size=10, face="bold"),
          axis.title.x=element_text(size=12, face="bold"),
          axis.title.y=element_text(size=12, face="bold"))+
    theme(strip.text.x = element_text(size=10, face="bold"),
          strip.text.y = element_text(size=10, face="bold"))+
    geom_vline(aes(xintercept = AllInsLenMean), colour="blue", size=0.75,linetype = "longdash")
  #geom_vline(aes(xintercept = AllInsMedian), colour="red", size=0.75,linetype = "longdash")
  AllInsertionLengthPlot
  pdf(paste(outdir, "/", "plots/", plasmid, "_Insertion_Length_Plot_All.pdf", sep=""), width = 20)
  AllInsertionLengthPlot
  dev.off()
  
  #--------------------------Insertion Length Plot - Consistent----------------
  ConsistentInsertion = subset(all, CLASS == "Insertion" & CONSISTENCY =="TRUE")
  ConsistentReadsAggInsertionLength = aggregate(READS~INSERTION_LENGTH, data=ConsistentInsertion, sum)
  ConsistentReadsAggInsertionLength$Percent_Insertion = ConsistentReadsAggInsertionLength$READS / sum(ConsistentReadsAggInsertionLength$READS)*100
  ConsistentReadsAggInsertionLength$READSxLength = ConsistentReadsAggInsertionLength$READS * ConsistentReadsAggInsertionLength$INSERTION_LENGTH
  ConsInsLenMean = sum(ConsistentReadsAggInsertionLength$READSxLength)/sum(ConsistentReadsAggInsertionLength$READS)
  
  ConsistentInsertionLengthPlot = ggplot(ConsistentReadsAggInsertionLength) + geom_bar(aes(x=INSERTION_LENGTH, y=Percent_Insertion), stat="identity", colour="black", fill="grey50")+
    theme_bw()+scale_y_continuous(name = "Percent of SD-MMEJ Consistent Insertions")+ scale_x_continuous(name = "Insertion Length (bp)", breaks = seq(1, max(ReadsAggInsertionLength$INSERTION_LENGTH), by = 1))+
    theme(plot.title = element_text(color="grey28", size=12),
          axis.text.x=element_text(size=10, face="bold"),
          axis.text.y=element_text(size=10, face="bold"),
          axis.title.x=element_text(size=12, face="bold"),
          axis.title.y=element_text(size=12, face="bold"))+
    theme(strip.text.x = element_text(size=10, face="bold"),
          strip.text.y = element_text(size=10, face="bold"))+
    geom_vline(aes(xintercept = ConsInsLenMean), colour="blue", size=0.75,linetype = "longdash")
  ConsistentInsertionLengthPlot
  pdf(paste(outdir, "/", "plots/", plasmid, "_Insertion_Length_Plot_Consistent.pdf", sep=""), width = 15)
  ConsistentInsertionLengthPlot
  dev.off()
}

#----------------Deletion Repeat Motif Plot - Data Manipulation-------------------
if (nrow(del) > 0){
  SampleID = as.data.frame(cbind(del[,2], del[,5]))
  colnames(SampleID)=c("Sample.ID", "RECONSTRUCTED_SEQ")
  
  ReconSeq = as.data.frame(sapply(SampleID$RECONSTRUCTED_SEQ, gsub, pattern="-", replacement=""))
  colnames(ReconSeq)[1]="RECONSTRUCTED_SEQ"
  SampleID$RECONSTRUCTED_SEQ=ReconSeq$RECONSTRUCTED_SEQ
  
  DeletionData2 = merge(x= DeletionData, y= SampleID, by.x="Sample.ID", by.y="Sample.ID")
  DeletionData2 = merge(DeletionData2, CombinedCurated, by="RECONSTRUCTED_SEQ", all.x=TRUE)
  DeletionData2$start = ifelse(DeletionData2$Break.Side=="left", (BreakPointFromLeft+1)-DeletionData2$Motif.to.Break,
                               (BreakPointFromLeft+1)+DeletionData2$Motif.to.Break-DeletionData2$Motif.Length)
  DeletionData2$end = ifelse(DeletionData2$Break.Side=="left", BreakPointFromLeft-DeletionData2$Motif.to.Break+DeletionData2$Motif.Length,
                             (BreakPointFromLeft+DeletionData2$Motif.to.Break))
  DeletionData2$motif_mechID = paste(DeletionData2$start, DeletionData2$end, DeletionData2$Mechanism, sep="-")
  DelAgg = aggregate(READS~motif_mechID, data=DeletionData2, sum)
  DelAggTab = as.data.frame(table(DeletionData2$motif_mechID))
  DelAgg = merge(DelAgg, DelAggTab, by.x="motif_mechID", by.y="Var1")
  DeletionData2 = merge(DelAgg, DeletionData2, by="motif_mechID")
  DeletionData2 = DeletionData2[!duplicated(DeletionData2$motif_mechID),]
  DeletionLeft = subset(DeletionData2, start<BreakPointFromLeft)
  DeletionRight = subset(DeletionData2, start>BreakPointFromLeft)
  DeletionLeft = DeletionLeft[order(DeletionLeft$start),]
  DeletionRight = DeletionRight[order(DeletionRight$start),]
  DeletionLeft$order = paste(10:(9+nrow(DeletionLeft)), "-", sep = "")
  DeletionRight$order = paste(10:(9+nrow(DeletionRight)), "-", sep="")
  DeletionData2 = rbind(DeletionLeft, DeletionRight)
  DeletionData2 = DeletionData2[order(DeletionData2$order, decreasing=TRUE),]
  DeletionData3 = DeletionData2$motif_mechID
  DeletionData2 = arrange(transform(DeletionData2, motif_mechID=factor(motif_mechID, levels=motif_mechID)),motif_mechID)
  DeletionData2$percent_deletion = DeletionData2$READS.x/sum(DeletionData2$READS.x)*100
  DeletionData2$MOTIF_START2 = as.numeric(DeletionData2$start)-0.5
  DeletionData2$MOTIF_END2 = as.numeric(DeletionData2$end)+0.5
  write.csv(DeletionData2, paste(outdir, "/", "table_outputs/", plasmid, "_del_data_for_template_plot.csv", sep=""))
  
  DelMotif = data.frame(x=numeric(), y=numeric(), z=numeric(), w=numeric(), stringsAsFactors = FALSE)
  for(i in 1:nrow(DeletionData2)){
    g = as.data.frame(DeletionData2[i, 48]:DeletionData2[i,49])
    DelMotif1 = cbind(DeletionData2[i,1],g,DeletionData2[i,2], DeletionData2[i, 8])
    DelMotif = rbind(DelMotif, DelMotif1)
  }
  colnames(DelMotif) = c("motif_mechID", "temp_coord","READS","mechanism")
  DelSnapBack = subset(DelMotif, mechanism=="snap-back")
  DelLoopOut = subset(DelMotif, mechanism=="loop-out")
  DelSnapBack = aggregate(READS~temp_coord, data=DelSnapBack, sum)
  DelSnapBack$mechanism = "Snap-back"
  DelLoopOut = aggregate(READS~temp_coord, data=DelLoopOut, sum)
  DelLoopOut$mechanism = "Loop-out"
  DelMotif2 = rbind(DelSnapBack, DelLoopOut)
  DelMotif2$x = "1"
  DelMotif2$percent_deletion = DelMotif2$READS/sum(DelMotif2$READS)*100
  DelMotif2$temp_coord2 = DelMotif2$temp_coord-0.5
  write.csv(DelMotif2, paste(outdir, "/", "table_outputs/", plasmid, "_del_data_for_temp_plot_mech.csv", sep=""))
  
  DelMotif = aggregate(READS~temp_coord, data=DelMotif, sum)
  DelMotif$x = plasmid
  DelMotif$percent_deletion = DelMotif$READS/sum(DelMotif$READS)*100
  DelMotif$temp_coord2 = DelMotif$temp_coord-0.5
  write.csv(DelMotif, paste(outdir, "/", "table_outputs/", plasmid, "_del_data_for_temp_plot.csv", sep=""))
  
  #----------------Deletion Repeat Motif Plot-------------------
  DeletionData2 = arrange(transform(DeletionData2, motif_mechID=factor(motif_mechID, levels = motif_mechID)), motif_mechID)
  
  DelMotifPlotBreaks = min(DeletionData2$MOTIF_START2-1.5):max(DeletionData2$MOTIF_END2+2.5)
  DelMotifPlotLabel = dput(as.character(referenceSplit[DelMotifPlotBreaks]))
  rect_left <- c(33.5,43.5,53.5,63.5,73.5,83.5,93.5,103.5,113.5,123.5,133.5,143.5,153.5,163.5,173.5,183.5,193.5,203.5,213.5,223.5)
  rectangles <- data.frame(
    xmin = rect_left,
    xmax = rect_left + 5,
    ymin = -Inf,
    ymax = Inf
  )
  DelMotifPlot = ggplot()+
    geom_boxplot(data=DeletionData2, aes(x=motif_mechID, ymin = MOTIF_START2, lower = MOTIF_START2, middle = MOTIF_START2, 
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
                         limits=c(0,max(DeletionData2$percent_deletion+1.5)),
                         guide_legend(title="% SD-MMEJ\nConsistent\nDeletion Reads"))+
    scale_y_continuous(limits=c(min(DeletionData2$MOTIF_START2-1.5),max(DeletionData2$MOTIF_END2+1.5)),
                       breaks=DelMotifPlotBreaks,
                       labels = DelMotifPlotLabel,
                       name = "Nucleotide",
                       expand=c(0,0))+
    theme_classic(base_size = 10)+
    theme(legend.text=element_text(size=8,face="bold"),
          axis.text.y= element_blank(),
          axis.title.y=element_blank(),
          axis.text.x= element_text(face="bold"),
          axis.title.x=element_text(face="bold"),
          axis.ticks=element_blank(),
          legend.key.size=unit(0.42,"cm"),
          panel.border = element_rect(colour = "black", fill=NA, size=1))+
    guides(guide_legend(title.theme = element_text(size=4, face="bold", angle=0),
                        label.theme=element_text(size=8, face="bold", angle=0)))+
    geom_rect(data=rectangles, aes(ymin=as.numeric(xmin), ymax=as.numeric(xmax), xmin=ymin, xmax=ymax),
              fill="grey", alpha=0.8, colour=NA)+
    geom_boxplot(data=DeletionData2, aes(x=motif_mechID, ymin = MOTIF_START2, lower = MOTIF_START2, middle = MOTIF_START2, 
                                         upper = MOTIF_END2, ymax = MOTIF_END2,
                                         fill=percent_deletion, linetype=Mechanism),
                 colour="black",
                 stat = "identity",
                 width=.8, 
                 lwd=.5)+  
    geom_hline(yintercept = BreakPointFromLeft+0.5, colour="red", size=0.75)
  DelMotifPlot
  pdf(paste(outdir, "/", "plots/", plasmid, "_Deletion_Repeat_Motif.pdf", sep=""), height = 7, width = 10)
  DelMotifPlot
  dev.off()
  
  #----------------Deletion Flap Plot - Data Manipulation-----------------------------------
  DelFlap = DeletionData2[,-c(1,3,9:16, 18, 20:27, 36, 39:43, 46,47)]
  DelFlap = subset(DelFlap, TOTAL_DELETION !="NA")
  DelFlap$LeftFlap = BreakPointFromLeft+0.5+DelFlap$DELETION_FROM_LEFT
  DelFlap$RightFlap = BreakPointFromLeft+0.5-DelFlap$DELETION_FROM_RIGHT
  DelFlap$FlapID = paste(DelFlap$LeftFlap, DelFlap$RightFlap, DelFlap$Mechanism, DelFlap$Repair.Type)
  DelFlapAggr=aggregate(READS.x~FlapID, data=DelFlap,sum)
  DelFlapTable=as.data.frame(table(DelFlap$FlapID))
  DelFlapAggr=merge(DelFlapAggr, DelFlapTable, by.x="FlapID", by.y="Var1")
  DelFlap=merge(DelFlapAggr, DelFlap, by="FlapID")
  DelFlap=DelFlap[!duplicated(DelFlap$FlapID),]
  colnames(DelFlap)[2]="READS"
  DelFlapMHJ = subset(DelFlap, Repair.Type=="MHJ")
  DelFlapABJ = subset(DelFlap, Repair.Type=="ABJ")
  DelFlapMHJ$percent_flap = DelFlapMHJ$READS/(sum(DelFlapMHJ$READS))*100
  DelFlapABJ$percent_flap = DelFlapABJ$READS/(sum(DelFlapABJ$READS))*100
  
  #----------------Deletion Flap Plot - MHJ Plot-----------------------------------
  rect_left <- c(43.5,53.5,63.5,73.5,83.5,93.5,103.5,113.5,123.5,133.5,143.5,153.5,163.5,173.5,183.5,193.5,203.5)
  DelFlapMHJPlotBreaks = min(DelFlapMHJ$LeftFlap-1.5):max(DelFlapMHJ$RightFlap+1.5)
  DelFlapMHJPlotLabel = dput(as.character(referenceSplit[DelFlapMHJPlotBreaks]))
  DelFlapMHJPlot <- ggplot()+
    geom_boxplot(data=DelFlapMHJ, aes(x=FlapID, ymin = LeftFlap, lower = LeftFlap, middle = LeftFlap, 
                                      upper = RightFlap, ymax = RightFlap,fill=percent_flap, linetype=Mechanism),
                 colour="white", stat = "identity",width=.8,lwd=.5) + coord_flip()+scale_linetype_manual(values=c("solid", "dashed","dotted"),guide="none")+
    scale_fill_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),values=c(1.0,0.8,0.6,0.4,0.2,0),
                         guide_legend(title="% Deletion From\nSingle Step\nMHJ Deletion Events"),limits=c(0,max(DelFlapMHJ$percent_flap+1.5)))+scale_y_continuous(limits=c(min(DelFlapMHJ$LeftFlap-1.5),max(DelFlapMHJ$RightFlap+1.5)),
                                                                                                                                                                  breaks=DelFlapMHJPlotBreaks,
                                                                                                                                                                  labels = DelFlapMHJPlotLabel,
                                                                                                                                                                  name = "Nucleotide",expand=c(0,0))+theme_classic(base_size = 10)+theme(legend.text=element_text(size=8,face="bold"),
                                                                                                                                                                                                                                         axis.text.y= element_blank(),axis.title.y=element_blank(),axis.text.x= element_text(face="bold"),axis.title.x=element_text(face="bold"),
                                                                                                                                                                                                                                         axis.ticks=element_blank(),
                                                                                                                                                                                                                                         panel.border = element_rect(colour = "black", fill=NA, size=1))+guides(guide_legend(title.theme = element_text(size=8, face="bold", angle=0),
                                                                                                                                                                                                                                                                                                                             label.theme=element_text(size=8, face="bold", angle=0)))+geom_rect(data=rectangles, aes(ymin=as.numeric(xmin), 
                                                                                                                                                                                                                                                                                                                                                                                                                     ymax=as.numeric(xmax), xmin=ymin, xmax=ymax),fill="grey", alpha=0.4, colour=NA)+
    geom_boxplot(data=DelFlapMHJ, aes(x=FlapID, ymin = LeftFlap, lower = LeftFlap, middle = LeftFlap, 
                                      upper = RightFlap, ymax = RightFlap,fill=percent_flap, linetype=Mechanism), colour="black",stat = "identity",
                 width=.8,lwd=.5)+ geom_hline(yintercept = BreakPointFromLeft+0.5, colour="red", size=0.75)
  
  DelFlapMHJPlot
  pdf(paste(outdir, "/", "plots/", plasmid, "_Deletion_Flap_MHJ_Plot.pdf", sep=""), height = 5, width = 10)
  DelFlapMHJPlot
  dev.off()
  
  #----------------Deletion Flap Plot - ABJ Plot-----------------------------------
  rect_left <- c(43.5,53.5,63.5,73.5,83.5,93.5,103.5,113.5,123.5,133.5,143.5,153.5,163.5,173.5,183.5,193.5,203.5)
  DelFlapABJPlotBreaks = min(DelFlapABJ$LeftFlap-1.5):max(DelFlapABJ$RightFlap+1.5)
  DelFlapABJPlotLabel = dput(as.character(referenceSplit[DelFlapABJPlotBreaks]))
  DelFlapABJPlot <- ggplot()+
    geom_boxplot(data=DelFlapABJ, aes(x=FlapID, ymin = LeftFlap, lower = LeftFlap, middle = LeftFlap, 
                                      upper = RightFlap, ymax = RightFlap,fill=percent_flap, linetype=Mechanism),
                 colour="white", stat = "identity",width=.8,lwd=.5) + coord_flip()+scale_linetype_manual(values=c("solid", "dashed","dotted"),guide="none")+
    scale_fill_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),values=c(1.0,0.8,0.6,0.4,0.2,0),
                         guide_legend(title="% Deletion From\nSingle Step\nABJ Deletion Events"),limits=c(0,max(DelFlapABJ$percent_flap+1.5)))+scale_y_continuous(limits=c(min(DelFlapABJ$LeftFlap-1.5),max(DelFlapABJ$RightFlap+1.5)),
                                                                                                                                                                  breaks=DelFlapABJPlotBreaks,
                                                                                                                                                                  labels = DelFlapABJPlotLabel,
                                                                                                                                                                  name = "Nucleotide",expand=c(0,0))+theme_classic(base_size = 10)+theme(legend.text=element_text(size=8,face="bold"),
                                                                                                                                                                                                                                         axis.text.y= element_blank(),axis.title.y=element_blank(),axis.text.x= element_text(face="bold"),axis.title.x=element_text(face="bold"),
                                                                                                                                                                                                                                         axis.ticks=element_blank(),
                                                                                                                                                                                                                                         panel.border = element_rect(colour = "black", fill=NA, size=1))+guides(guide_legend(title.theme = element_text(size=8, face="bold", angle=0),
                                                                                                                                                                                                                                                                                                                             label.theme=element_text(size=8, face="bold", angle=0)))+geom_rect(data=rectangles, aes(ymin=as.numeric(xmin), 
                                                                                                                                                                                                                                                                                                                                                                                                                     ymax=as.numeric(xmax), xmin=ymin, xmax=ymax),fill="grey", alpha=0.4, colour=NA)+
    geom_boxplot(data=DelFlapABJ, aes(x=FlapID, ymin = LeftFlap, lower = LeftFlap, middle = LeftFlap, 
                                      upper = RightFlap, ymax = RightFlap,fill=percent_flap, linetype=Mechanism), colour="black",stat = "identity",
                 width=.8,lwd=.5)+ geom_hline(yintercept = BreakPointFromLeft+0.5, colour="red", size=0.75)
  
  DelFlapABJPlot
  pdf(paste(outdir, "/", "plots/", plasmid, "_Deletion_Flap_ABJ_Plot.pdf", sep=""), height = 5, width = 10)
  DelFlapABJPlot
  dev.off()
  
  #----------------Deletion Resection - Data Manipulation-----------
  #To generate a plot of the furthest identifiable point of resection or duplex unwinding = furthest from break repeat motif
  DelLeft=subset(DeletionData2, Break.Side == "left")
  DelLeft$RM_to_break= abs(DelLeft$P1.to.Break)
  DelRight=subset(DeletionData2, Break.Side=="right")
  DelRight$RM_to_break = abs(DelRight$P1.to.Break)
  DelFromBreak=rbind(DelLeft, DelRight)
  DelFromBreak=DelFromBreak[,-c(1,3, 20:23, 39, 40, 42,43,46, 47)]
  DelFromBreak$percent = (DelFromBreak$READS.x/sum(DelFromBreak$READS.x))*100
  write.csv(DelFromBreak, paste(outdir, "/", "table_outputs/", plasmid, "_deletion_resection_data.csv", sep=""))
  DelFromBreakMean=sum(DelFromBreak$RM_to_break*DelFromBreak$READS.x)/sum(DelFromBreak$READS.x)
  
  #----------------Deletion Resection - Plot - Break Side-----------
  DeletionResectionPlotSide = ggplot(DelFromBreak)+
    geom_bar(aes(x=RM_to_break, y=percent, color=Break.Side), position="stack", stat="identity", fill="grey50", width=0.5) +
    theme_bw()+
    geom_vline(aes(xintercept = DelFromBreakMean), colour="blue", size=0.75,linetype = "longdash")+
    #geom_vline(aes(xintercept = median), colour="red", size=0.75,linetype = "longdash",data=distMedian)+
    scale_y_continuous(name = "Percent Deletions (%)")+  
    scale_x_continuous(name="Resection/Unwinding Required for SD-MMEJ Action (bp)",
                       breaks = seq(0, max(DelFromBreak$RM_to_break), by = 2))+
    theme(axis.text.x=element_text(size=10, face="bold"),
          axis.text.y=element_text(size=10, face="bold"),
          axis.title.x=element_text(size=12, face="bold"),
          axis.title.y=element_text(size=12, face="bold"))+
    theme(strip.text.x = element_text(size=10, face="bold"),
          strip.text.y = element_text(size=10, face="bold"))
  DeletionResectionPlotSide
  pdf(paste(outdir, "/", "plots/", plasmid, "_Deletion_Resection_Plot_BreakSide.pdf", sep=""), height = 5, width = 10)
  DeletionResectionPlotSide
  dev.off()
  
  #----------------Deletion Resection - Plot - Mechanism-----------
  DeletionResectionPlotMechanism = ggplot(DelFromBreak)+
    geom_bar(aes(x=RM_to_break, y=percent, color=Mechanism), position="stack", stat="identity", fill="grey50", width=0.5) +
    theme_bw()+
    geom_vline(aes(xintercept = DelFromBreakMean), colour="blue", size=0.75,linetype = "longdash")+
    #geom_vline(aes(xintercept = median), colour="red", size=0.75,linetype = "longdash",data=distMedian)+
    scale_y_continuous(name = "Percent Deletions (%)")+  
    scale_x_continuous(name="Resection/Unwinding Required for SD-MMEJ Action (bp)",
                       breaks = seq(0, max(DelFromBreak$RM_to_break), by = 2))+
    theme(axis.text.x=element_text(size=10, face="bold"),
          axis.text.y=element_text(size=10, face="bold"),
          axis.title.x=element_text(size=12, face="bold"),
          axis.title.y=element_text(size=12, face="bold"))+
    theme(strip.text.x = element_text(size=10, face="bold"),
          strip.text.y = element_text(size=10, face="bold"))
  DeletionResectionPlotMechanism
  pdf(paste(outdir, "/", "plots/", plasmid, "_Deletion_Resection_Plot_Mechanism.pdf", sep=""), height = 5, width = 10)
  DeletionResectionPlotMechanism
  dev.off()
  
  #----------------Deletion Resection - Plot - Repair Type-----------
  DeletionResectionPlotRepairType = ggplot(DelFromBreak)+
    geom_bar(aes(x=RM_to_break, y=percent, color=Repair.Type), position="stack", stat="identity", fill="grey50", width=0.5) +
    theme_bw()+
    geom_vline(aes(xintercept = DelFromBreakMean), colour="blue", size=0.75,linetype = "longdash")+
    #geom_vline(aes(xintercept = median), colour="red", size=0.75,linetype = "longdash",data=distMedian)+
    scale_y_continuous(name = "Percent Deletions (%)")+  
    scale_x_continuous(name="Resection/Unwinding Required for SD-MMEJ Action (bp)",
                       breaks = seq(0, max(DelFromBreak$RM_to_break), by = 2))+
    theme(axis.text.x=element_text(size=10, face="bold"),
          axis.text.y=element_text(size=10, face="bold"),
          axis.title.x=element_text(size=12, face="bold"),
          axis.title.y=element_text(size=12, face="bold"))+
    theme(strip.text.x = element_text(size=10, face="bold"),
          strip.text.y = element_text(size=10, face="bold"))
  DeletionResectionPlotRepairType
  pdf(paste(outdir, "/", "plots/", plasmid, "_Deletion_Resection_Plot_Repair_Type.pdf", sep=""), height = 5, width = 10)
  DeletionResectionPlotRepairType
  dev.off()
  
  #----------------Deletion Repeat Motif by Side - Plot - Mechanism---------------
  DelSide = DeletionData2
  DelSideMechPlot = ggplot(DelSide)+
    geom_bar(aes(x=Break.Side, y=percent_deletion, color=Mechanism), position="stack", stat="identity", fill="grey50") +
    theme_bw()+
    scale_y_continuous(name = "Percent Usage (%)")+  
    scale_x_discrete(name="Side")+
    theme(axis.text.x=element_text(size=10, face="bold"),
          axis.text.y=element_text(size=10, face="bold"),
          axis.title.x=element_text(size=12, face="bold"),
          axis.title.y=element_text(size=12, face="bold"))+
    theme(strip.text.x = element_text(size=10, face="bold"),
          strip.text.y = element_text(size=10, face="bold"))
  DelSideMechPlot
  pdf(paste(outdir, "/", "plots/", plasmid, "_Deletion_Side_Usage_Mech_Plot.pdf", sep=""), height = 5, width = 10)
  DelSideMechPlot
  dev.off()
  
  #----------------Deletion Repeat Motif by Side - Plot - Repair Type---------------
  DelSideTypePlot = ggplot(DelSide)+
    geom_bar(aes(x=Break.Side, y=percent_deletion, color=Repair.Type), position="stack", stat="identity", fill="grey50") +
    theme_bw()+
    scale_y_continuous(name = "Percent Usage (%)")+  
    scale_x_discrete(name="Side")+
    theme(axis.text.x=element_text(size=10, face="bold"),
          axis.text.y=element_text(size=10, face="bold"),
          axis.title.x=element_text(size=12, face="bold"),
          axis.title.y=element_text(size=12, face="bold"))+
    theme(strip.text.x = element_text(size=10, face="bold"),
          strip.text.y = element_text(size=10, face="bold"))
  DelSideTypePlot
  pdf(paste(outdir, "/", "plots/", plasmid, "_Deletion_Side_Usage_Repair_Type_Plot.pdf", sep=""), height = 5, width = 10)
  DelSideTypePlot
  dev.off()
  
  
  ## second set of deletion plots
  
  #----------------------------Microhomology Length Plot - Inaccurate Reads----------
  MHevents = subset(all, all$MH_Length>0)
  MHlength = aggregate(as.numeric(percent_inaccurate)~MH_Length, data=MHevents, sum)
  colnames(MHlength)[2]="percent_inaccurate"
  
  MHLengthPlot = ggplot(MHlength)+
    geom_bar(aes(x=MH_Length, y=percent_inaccurate), position="stack", stat="identity", colour="black", fill="grey50") +
    theme_bw()+
    #geom_vline(aes(xintercept = mean), colour="blue", size=0.75,linetype = "longdash",data=len.mean)+
    #geom_vline(aes(xintercept = median), colour="red", size=0.75,linetype = "longdash",data=len.median)+
    scale_y_continuous(name = "Percent Inaccurate Reads")+
    scale_x_continuous(name="Microhomology (bp)",
                       breaks = seq(1, 100, by = 1))+
    theme(axis.text.x=element_text(size=10, face="bold"),
          axis.text.y=element_text(size=10, face="bold"),
          axis.title.x=element_text(size=12, face="bold"),
          axis.title.y=element_text(size=12, face="bold"))+
    theme(strip.text.x = element_text(size=10, face="bold"),
          strip.text.y = element_text(size=10, face="bold"))
  MHLengthPlot
  pdf(paste(outdir, "/", "plots/", plasmid, "_MH_Length_Plot_Consistent_InaccurateReads.pdf", sep=""), width = 10)
  MHLengthPlot
  dev.off()
  
  #----------------------------Microhomology Length Plot - Deletion Events----------
  MHevents = subset(all, all$MH_Length>0)
  MHevents$percent_deletion = MHevents$READS/sum(MHevents$READS)*100
  MHlengthdeletion = aggregate(as.numeric(percent_deletion)~MH_Length, data=MHevents, sum)
  colnames(MHlengthdeletion)[2]="percent_deletion"
  
  MHLengthPlotDeletion = ggplot(MHlengthdeletion)+
    geom_bar(aes(x=MH_Length, y=percent_deletion), position="stack", stat="identity", colour="black", fill="grey50") +
    theme_bw()+
    #geom_vline(aes(xintercept = mean), colour="blue", size=0.75,linetype = "longdash",data=len.mean)+
    #geom_vline(aes(xintercept = median), colour="red", size=0.75,linetype = "longdash",data=len.median)+
    scale_y_continuous(name = "Percent SD-MMEJ Consistent Deletions")+
    scale_x_continuous(name="Microhomology (bp)",
                       breaks = seq(1, 100, by = 1))+
    theme(axis.text.x=element_text(size=10, face="bold"),
          axis.text.y=element_text(size=10, face="bold"),
          axis.title.x=element_text(size=12, face="bold"),
          axis.title.y=element_text(size=12, face="bold"))+
    theme(strip.text.x = element_text(size=10, face="bold"),
          strip.text.y = element_text(size=10, face="bold"))
  MHLengthPlotDeletion
  pdf(paste(outdir, "/", "plots/", plasmid, "_MH_Length_Plot_Deletion_Consistent.pdf", sep=""), width = 10)
  MHLengthPlot
  dev.off()
  
  #----------------------------Microhomology Usage Plot- Inaccurate Reads----------
  MHusage = aggregate(as.numeric(percent_inaccurate)~MICROHOMOLOGY, data=MHevents, sum)
  colnames(MHusage)[2]="percent_inaccurate"
  microhomology = MHusage$MICROHOMOLOGY
  
  MHusagePlot <- ggplot(MHusage)+
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
  MHusagePlot
  pdf(paste(outdir, "/", "plots/", plasmid, "_MH_Usage_Plot_InaccurateReads.pdf", sep=""), width = 15)
  MHusagePlot
  dev.off()
  
  #----------------------------Microhomology Usage Plot - Deletion Events----------
  MHusageDeletion = aggregate(as.numeric(percent_deletion)~MICROHOMOLOGY, data=MHevents, sum)
  colnames(MHusageDeletion)[2]="percent_deletion"
  microhomology = MHusageDeletion$MICROHOMOLOGY
  write.csv(MHusageDeletion, paste(outdir, "/", "table_outputs/", plasmid, "_MH_usage_consistent_deletions.csv", sep=""))
  
  MHusagePlotDeletion <- ggplot(MHusageDeletion)+
    geom_histogram(aes(x=MICROHOMOLOGY, y=percent_deletion), position="stack", stat="identity", colour="black", fill="grey50") +
    theme_bw()+
    # facet_wrap(~PG, scales="free_y")+
    # geom_vline(aes(xintercept = mean), colour="blue", size=0.75,linetype = "longdash",data=len.mean)+
    # geom_vline(aes(xintercept = median), colour="red", size=0.75,linetype = "longdash",data=len.median)+
    # ggtitle("Left Synthesis")+
    scale_y_continuous(name = "Percent SD-MMEJ Consistent Deletions")+  
    scale_x_discrete(name="Microhomology")+
    theme(plot.title = element_text(color="grey28", size=12),
          axis.text.x=element_text(size=10, face="bold"),
          axis.text.y=element_text(size=10, face="bold"),
          axis.title.x=element_text(size=12, face="bold"),
          axis.title.y=element_text(size=12, face="bold"))+
    theme(strip.text.x = element_text(size=10, face="bold"),
          strip.text.y = element_text(size=10, face="bold"))
  MHusagePlotDeletion
  pdf(paste(outdir, "/", "plots/", plasmid, "_MH_Usage_Plot_Consistent_Deletions.pdf", sep=""), width = 15)
  MHusagePlotDeletion
  dev.off()
  
  
  #--------------------------Deletion Length Plot - All----------------
  deletions$deletion_length = abs(deletions$TOTAL_DELETION)
  AllDeletionsAggLength = aggregate(READS~deletion_length + CONSISTENCY, data=deletions, sum)
  AllDeletionsAggLength$Percent_Deletion = AllDeletionsAggLength$READS/sum(AllDeletionsAggLength$READS)*100
  AllDeletionsAggLength$READSxLength = AllDeletionsAggLength$READS * AllDeletionsAggLength$deletion_length
  AllDeletionsAggMean = sum(AllDeletionsAggLength$READSxLength)/sum(AllDeletionsAggLength$READS)
  
  AllDeletionsLengthPlot = ggplot(AllDeletionsAggLength)+geom_bar(aes(x=deletion_length, y=Percent_Deletion, fill= CONSISTENCY), stat="identity", colour="black")+
    theme_bw()+scale_y_continuous(name = "Percent of Deletions")+ scale_x_continuous(name = "Deletion Length (bp)", breaks = seq(0, max(AllDeletionsAggLength$deletion_length), by =2))+
    theme(plot.title = element_text(color="grey28", size=12),
          axis.text.x=element_text(size=10, face="bold"),
          axis.text.y=element_text(size=10, face="bold"),
          axis.title.x=element_text(size=12, face="bold"),
          axis.title.y=element_text(size=12, face="bold"))+
    theme(strip.text.x = element_text(size=10, face="bold"),
          strip.text.y = element_text(size=10, face="bold"))+
    geom_vline(aes(xintercept = AllDeletionsAggMean), colour="blue", size=0.75,linetype = "longdash")
  AllDeletionsLengthPlot
  pdf(paste(outdir, "/", "plots/", plasmid, "_Deletion_Length_Plot_All.pdf", sep=""), width = 20)
  AllDeletionsLengthPlot
  dev.off()
  
  #--------------------------Deletion Length Plot - Consistent----------------
  ConsistentDeletions = subset(deletions, CONSISTENCY=="TRUE")
  ConsistentDeletionsAggLength = aggregate(READS~deletion_length, data=ConsistentDeletions, sum)
  ConsistentDeletionsAggLength$Percent_Deletion = ConsistentDeletionsAggLength$READS/sum(ConsistentDeletionsAggLength$READS)*100
  ConsistentDeletionsAggLength$READSxLength = ConsistentDeletionsAggLength$READS * ConsistentDeletionsAggLength$deletion_length
  ConDelAggMean = sum(ConsistentDeletionsAggLength$READSxLength)/sum(ConsistentDeletionsAggLength$READS)
  
  ConsistentDeletionsLengthPlot = ggplot(ConsistentDeletionsAggLength)+geom_bar(aes(x=deletion_length, y=Percent_Deletion), stat="identity", colour="black", fill="grey50")+
    theme_bw()+scale_y_continuous(name = "Percent of SD-MMEJ Consistent Deletions")+ scale_x_continuous(name = "Deletion Length (bp)", breaks = seq(0, max(ConsistentDeletionsAggLength$deletion_length), by = 2))+
    theme(plot.title = element_text(color="grey28", size=12),
          axis.text.x=element_text(size=10, face="bold"),
          axis.text.y=element_text(size=10, face="bold"),
          axis.title.x=element_text(size=12, face="bold"),
          axis.title.y=element_text(size=12, face="bold"))+
    theme(strip.text.x = element_text(size=10, face="bold"),
          strip.text.y = element_text(size=10, face="bold"))+
    geom_vline(aes(xintercept = ConDelAggMean), colour="blue", size=0.75,linetype = "longdash")
  ConsistentDeletionsLengthPlot
  pdf(paste(outdir, "/", "plots/", plasmid, "_Deletion_Length_Plot_Consistent.pdf", sep=""), width = 20)
  ConsistentDeletionsLengthPlot
  dev.off()
}

#----------------Consistency and Inaccurate Repair Plots - Data Manipulation------------------------
#use the "all" Data Frame or import the _all_SD-MMEJ_consistency.csv
#all = read.csv("table_outputs/_all_SD-MMEJ_consistency.csv", sep = ",", header = TRUE, row.names=1)

## THIS IS THE SAME PROBLEM WITH "all" colnames
view(all)
Consistent = subset(all, CONSISTENCY=="TRUE")
head()
ConsistentReads = sum(Consistent$READS)
Inconsistent = subset(all, CONSISTENCY=="FALSE")
write.csv(Inconsistent, paste(outdir, "/", "table_outputs/", plasmid, "_Inconsistent_data.csv", sep=""))
      #Open the .csv just written and sort through the inconsistent large complex events for real versus sequencing artifacts
Inconsistent= read.csv(paste(outdir, "/", "table_outputs/", plasmid, "_Inconsistent_data.csv", sep=""), row.names = 1)
InconsistentReads = sum(Inconsistent$READS)
AllInaccurate= rbind(Consistent, Inconsistent)
head(AllInaccurate)
AllInaccurate$percent_inaccurate=AllInaccurate$READS/sum(AllInaccurate$READS)*100
Consistent=subset(AllInaccurate, CONSISTENCY=="TRUE")
Inconsistent = subset(AllInaccurate, CONSISTENCY=="FALSE")
RepairType = c("InDel", "MHJ", "ABJ")
PercentOfInaccurate = aggregate(as.numeric(percent_inaccurate)~REPAIR_TYPE+plasmid,data=AllInaccurate, sum) #Breaks down the types of repair (consistent or not) amongst inaccurate repair events

head(AllInaccurate)


colnames(PercentOfInaccurate)[3]="percent_inaccurate"
PercentInaccurate = aggregate(as.numeric(percent_inaccurate)~plasmid,data=Consistent, sum) #Breaks down percent inaccurate repair that is consistent
colnames(PercentInaccurate)[2]="percent_inaccurate"
Indel2 = subset(PercentOfInaccurate, REPAIR_TYPE=="InDel")
MHJ2 = subset(PercentOfInaccurate, REPAIR_TYPE=="MHJ")
ABJ2 = subset(PercentOfInaccurate, REPAIR_TYPE=="ABJ")
RepairTypeOrdered = rbind(Indel2, MHJ2, ABJ2) # ordering in InDel MHJ ABJ order to match all previous presentations of this data

## If there aren't insertions the InDels will be empty
InDel = subset(AllInaccurate, REPAIR_TYPE=="InDel")
InDelReads = sum(as.numeric(InDel$READS))
MHJ = subset(AllInaccurate, REPAIR_TYPE=="MHJ")
MHJReads = sum(as.numeric(MHJ$READS))
ABJ = subset(AllInaccurate, REPAIR_TYPE=="ABJ")
ABJReads = sum(ABJ$READS)

InDelTrue = subset(InDel, CONSISTENCY=="TRUE")
InDelTrueReads = sum(InDelTrue$READS)
InDelFalse = subset(InDel, CONSISTENCY=="FALSE")
InDelFalseReads = sum(InDelFalse$READS)
InDelPerConsistent = InDelTrueReads / InDelReads *100

MHJTrue = subset(MHJ, CONSISTENCY=="TRUE")
MHJTrueReads = sum(MHJTrue$READS)
MHJFalse = subset(MHJ, CONSISTENCY=="FALSE")
MHJFalseReads = sum(MHJFalse$READS)
MHJPerConsistent = MHJTrueReads / MHJReads *100

ABJTrue = subset(ABJ, CONSISTENCY=="TRUE")
ABJTrueReads = sum(ABJTrue$READS)
ABJFalse = subset(ABJ, CONSISTENCY=="FALSE")
ABJFalseReads = sum(ABJFalse$READS)
ABJPerConsistent = ABJTrueReads / ABJReads *100

ConsistentRepairType = NULL
ConsistentRepairType$plasmid = c(plasmid, plasmid, plasmid)
ConsistentRepairType = as.data.frame(ConsistentRepairType)
ConsistentRepairType$REPAIR_TYPE = c("InDel", "MHJ", "ABJ")
ConsistentRepairType$percent_consistent = c(InDelPerConsistent, MHJPerConsistent, ABJPerConsistent)


#----------------Consistency and Inaccurate Repair Plots------------------------
InaccurateRepairPlot = ggplot(PercentOfInaccurate, aes(x=plasmid, y=percent_inaccurate, fill=REPAIR_TYPE, label=REPAIR_TYPE))+ 
  geom_bar(position="stack", stat="identity", colour="black")+
  theme_classic(base_size = 9)+ geom_text(position = position_stack(vjust=0.5), color=c("white","black","black"))+
  scale_fill_grey(start=0.0, end=1) +guides(fill="none")+
  scale_y_continuous(name = "% Inaccurate Reads", labels=c("0","25","50","75","100"), expand=c(0,0),limits=c(0,105),breaks=c(0,25,50,75,100))+
  scale_x_discrete(name = "plasmid")+
  theme(axis.text.x=element_text(size=9, face="bold"),axis.text.y=element_text(size=9, face="bold"),axis.title.x=element_text(size=10, face="bold"),
  axis.title.y=element_text(size=10, face="bold"))
InaccurateRepairPlot
pdf(paste(outdir, "/", "plots/", plasmid, "_Inaccurate_Repair_Plot_All.pdf",sep=""), width=5)
InaccurateRepairPlot
dev.off()

ConsistencyPlot = ggplot(PercentInaccurate, aes(fill=plasmid))+ 
  geom_bar(aes(x=plasmid, y=percent_inaccurate), width=0.67,stat="identity",colour="black") +
  ylab("% SD-MMEJ consistent")+
  scale_y_continuous(labels=c("0","25","50","75","100"),expand=c(0,0),limits=c(0,101),breaks=c(0,25,50,75,100))+xlab("plasmid")+
  scale_fill_manual(values=c("black", "gray100", "gray75", "gray50"))+guides(fill="none")+theme_classic(base_size = 9)+
  theme(axis.text.x=element_text(size=9, face="bold"),axis.text.y=element_text(size=9, face="bold"),axis.title.x=element_text(size=10, face="bold"),
  axis.title.y=element_text(size=10, face="bold"))
ConsistencyPlot
pdf(paste(outdir, "/", "plots/", plasmid, "_SD-MMEJ_Consistency_Plot.pdf", sep=""), width=5)
ConsistencyPlot
dev.off()

ConsistencyBreakdownPlot =  ggplot(ConsistentRepairType, aes(fill=plasmid))+ 
  geom_bar(aes(x=REPAIR_TYPE, y=percent_consistent), width=0.67,position = position_dodge(width=0.75),stat="identity",colour="black") +
  scale_x_discrete(limits=RepairType,expand=c(0,0.5),labels=c("InDel", "MHJ", "ABJ")) +ylab("% SD-MMEJ consistent")+
  scale_y_continuous(labels=c("0","25","50","75","100"),expand=c(0,0),limits=c(0,105),breaks=c(0,25,50,75,100))+
  xlab("Repair Junction")+guides(fill="none")+scale_fill_manual(values=c("black", "gray100", "gray75", "gray50"))+
  theme_classic(base_size = 9)+theme(legend.position=c(0,1), legend.justification=c(.1,.8),legend.key.size=unit(.25, "cm"))+
  theme(axis.text.x=element_text(size=9, face="bold"),axis.text.y=element_text(size=9, face="bold"),axis.title.x=element_text(size=10, face="bold"),
  axis.title.y=element_text(size=10, face="bold"))
ConsistencyBreakdownPlot
pdf(paste(outdir, "/", "plots/", plasmid, "_SD-MMEJ_Consistency_Breakdown_Plot.pdf", sep=""))
ConsistencyBreakdownPlot
dev.off()




#-------------------------Deletion Stopper Plot - Data Manipulation----------
#Project for another day

