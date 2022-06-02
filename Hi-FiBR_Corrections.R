HiFiBR=read.csv("Iw7_output/Iw7_reclassified.csv")
InsertionsOuput=read.csv("Iw7_output/Iw7_insertion_insertion_consistency2.csv")
ComplexOutput=read.csv("Iw7_output/Iw7_complex_insertion_consistency2.csv")
DeletionOutput=read.csv("Iw7_output/Iw7_deletion_consistency_table.txt", sep="\t", header=T)
DeletionSubset=read.csv("Iw7_output/Iw7_deletion_consistency_log_subset.csv")
names(HiFiBR)[names(HiFiBR)=="ALIGNED_SEQ"] = "RECONSTRUCTED_SEQ"
HiFiBR=subset(HiFiBR, READS>9)
inaccurate=subset(HiFiBR, CLASS!="exact")
exact=subset(HiFiBR, CLASS=="exact")
Distance_From_Break_Left=exact$INSERTION_START
Distance_From_Break_Right=exact$INSERTION_END
Deletions=subset(HiFiBR, CLASS=="deletion")
Complex=subset(HiFiBR, CLASS=="complex")
Insertions=subset(HiFiBR, CLASS="insertion")
Insertions$INSERTION_START_TRUE=Insertions$INSERTION_START+1 #defines where inserted DNA begins, from the left
Insertions$INSERTION_END_TRUE=Insertions$INSERTION_START+Insertions$INSERTION_LENGTH #defines where inserted DNA ends, from the left
Insertions$DELETION_START_TRUE = "NA" #defines where deleted DNA in reference starts, from the left
Insertions$DELETION_END_TRUE = "NA" #defines where deleted DNA in reference ends, from left

Deletions$INSERTION_START_TRUE="NA" #defines where inserted DNA begins, from the left
Deletions$INSERTION_END_TRUE="NA" #defines where inserted DNA ends, from the left
Deletions$START_TRUE=Distance_From_Break_Left+Deletions$DELETION_FROM_LEFT+1 #defines where deleted DNA in reference begins, from the left
Deletions$END_TRUE=Distance_From_Break_Left-Deletions$DELETION_FROM_RIGHT #defines where deleted DNAin reference ends, from the left

Complex$DELETION_START_TRUE=Distance_From_Break_Left+Complex$DELETION_FROM_LEFT+1
Complex$DELETION_END_TRUE=Distance_From_Break_Left-Complex$DELETION_FROM_RIGHT
Complex$INSERTION_START_TRUE=Complex$INSERTION_START+1
Complex$INSERTION_END_TRUE=Complex$INSERTION_START+Complex$INSERTION_LENGTH


#-----------Ignore everything below for now-------------
#-----------Next chunk of code is from plots script data organization-----------
InsertionsAndComplex = rbind(Complex, Insertions)
InsComOutput=rbind(InsertionsOuput, ComplexOutput)
InsComMerge = merge(InsComOutput,InsertionsAndComplex, "RECONSTRUCTED_SEQ", all=FALSE)
InsComMerge = InsComMerge[,-c(2:7,11, 13:20, 38:39)]
InsComConsistent=subset(InsComMerge, consistency=="TRUE")
InsComInconsistent=subset(InsComMerge, consistency=="FALSE")
names(DeletionSubset)[names(DeletionSubset)=="SEQ"] = "RECONSTRUCTED_SEQ"
DeletionSubset$RECONSTRUCTED_SEQ = toupper(DeletionSubset$RECONSTRUCTED_SEQ)
DeletionsMerged = merge(DeletionSubset, Deletions, "RECONSTRUCTED_SEQ", all=FALSE)


#----------Trying to remove obvious sequencing artifacts
InsComInconsistentArtifacts=subset(InsComInconsistent, InsComInconsistent$INSERTION_LENGTH>14)
SequencingArtifactsRight=subset(InsComInconsistentArtifacts, InsComInconsistentArtifacts$MATCH_RIGHT>Distance_From_Break_Right)
SequencingArtifactsLeft=subset(InsComInconsistentArtifacts, InsComInconsistentArtifacts$MATCH_LEFT>Distance_From_Break_Left)


#-----------Everything below is my old code to use as reference---------------
SequencingArtifactsRight=subset(Complex, Complex$MATCH_RIGHT>Distance_From_Break_Right)
SequencingArtifactsRight=subset(SequencingArtifactsRight, SequencingArtifactsRight$MATCH_LEFT+SequencingArtifactsRight$MATCH_RIGHT==(Distance_From_Break_Left+Distance_From_Break_Right)-1)
SequencingArtifactsLeft=subset(Complex, Complex$MATCH_LEFT>Distance_From_Break_Left)
SequencingArtifactsLeft=subset(SequencingArtifactsLeft, SequencingArtifactsLeft$MATCH_LEFT+SequencingArtifactsLeft$MATCH_RIGHT==(Distance_From_Break_Left+Distance_From_Break_Right)-1)


