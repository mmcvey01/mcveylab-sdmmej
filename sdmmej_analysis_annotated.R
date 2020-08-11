# 07312017
# Varandt Y. Khodaverdian
#SD-MMEJ analysis script


# combined sam files from Steve Roberts into csv file. Replicates are not combined, so we will probably have to do that before we apply our cuttoff of 30 reads.
# experimental setup is 2 new mutations on the left side of the break, in S21 (mutant) and S22 (wt). These are the exp. We also injected into oregonR (con).
# So lets load it up!

rm(list=ls())
library("seqinr", lib.loc="~/R/win-library/3.1")
library(ggplot2)
library(grid)
library(gtable)
library(stringr)


#Working Directory
setwd("R:/Varandt/SD-MMEJ/07132017_amplicon_sequence")

a<-read.csv("combined.csv")
head(a)
tail(a)
m <- ggplot(a, aes(x = -TOTAL_DELETION, colour = SAMPLE))
m + geom_density(alpha=0.01, size=1) + geom_bar(stat='density', alpha=0.01,position="identity")+ facet_grid(PLASMID~ .)
# Looks cool
# Lets look at it with a >30 cuttof

# Need to combine replicates and add together the number of reads for matched replicates
a$RECONSTRUCTED_SEQ <-  gsub("-","",a$ALIGNED_SEQ) #remove dashes from aligned seq to make "reconstructed_seq" to use to find identical junctions
# seperate into appropriate groups to search
exp <- subset(a, SAMPLE=="exp")
con <- subset(a, SAMPLE=="con")
#seperate based on geneotype
expmut <- subset(exp, GENOTYPE=="mutant")
expwt <- subset(exp,GENOTYPE=="wt")
#seperate base on plasmid
expmutm4 <- subset(expmut,PLASMID=="M4")
expmutm5 <- subset(expmut,PLASMID=="M5")
expwtm4 <- subset(expwt,PLASMID=="M4")
expwtm5 <- subset(expwt,PLASMID=="M5")
conm4 <- subset(con,PLASMID=="M4")
conm5 <- subset(con,PLASMID=="M5")
#subset based on replicate so that we can merge then sum
expmutm4a <- subset(expmutm4,REP=="A")
expmutm4b <- subset(expmutm4,REP=="B")
expmutm5a <- subset(expmutm5,REP=="A")
expmutm5b <- subset(expmutm5,REP=="B")
expwtm4a <- subset(expwtm4,REP=="A")
expwtm4b <- subset(expwtm4,REP=="B")
expwtm5a <- subset(expwtm5,REP=="A")
expwtm5b <- subset(expwtm5,REP=="B")
conm4a <- subset(conm4, REP=="A")
conm4b <- subset(conm4, REP=="B")
conm5a <- subset(conm5, REP=="A")
conm5b <- subset(conm5, REP=="B")
#merge based on reconstructed seq
expmutm4merge <- merge(expmutm4a,expmutm4b, by="RECONSTRUCTED_SEQ",all=TRUE)
expmutm5merge <- merge(expmutm5a,expmutm5b, by="RECONSTRUCTED_SEQ",all=TRUE)
expwtm4merge <- merge(expwtm4a,expwtm4b, by="RECONSTRUCTED_SEQ",all=TRUE)
expwtm5merge <- merge(expwtm5a,expwtm5b, by="RECONSTRUCTED_SEQ",all=TRUE)
conm4merge <- merge(conm4a,conm4b, by="RECONSTRUCTED_SEQ",all=TRUE)
conm5merge <- merge(conm5a,conm5b, by="RECONSTRUCTED_SEQ",all=TRUE)
expmutm4merge[is.na(expmutm4merge)] <- 0
expmutm5merge[is.na(expmutm5merge)] <- 0
expwtm4merge[is.na(expwtm4merge)] <- 0
expwtm5merge[is.na(expwtm5merge)] <- 0
conm4merge[is.na(conm4merge)] <- 0
conm5merge[is.na(conm5merge)] <- 0
# sum the number of reads from A and B into total reads
expmutm4merge$total_reads <- expmutm4merge$X._OF_READS.x+expmutm4merge$X._OF_READS.y
expmutm5merge$total_reads <- expmutm5merge$X._OF_READS.x+expmutm5merge$X._OF_READS.y
expwtm4merge$total_reads <- expwtm4merge$X._OF_READS.x+expwtm4merge$X._OF_READS.y
expwtm5merge$total_reads <- expwtm5merge$X._OF_READS.x+expwtm5merge$X._OF_READS.y
conm4merge$total_reads <- conm4merge$X._OF_READS.x+conm4merge$X._OF_READS.y
conm5merge$total_reads <- conm5merge$X._OF_READS.x+conm5merge$X._OF_READS.y
# create mini table of RECONSTRUCTED_SEQ and total_reads to merge with sample tables to clean up
#m4 S21
mutm4.mini <- as.data.frame(expmutm4merge[,1])
mutm4.mini <- cbind(mutm4.mini,total_reads=expmutm4merge$total_reads)
names(mutm4.mini)[names(mutm4.mini) == 'expmutm4merge[, 1]'] <- 'RECONSTRUCTED_SEQ'
mutm4.mini <- unique(mutm4.mini)
expmutm4_test<- merge(expmutm4,mutm4.mini, by="RECONSTRUCTED_SEQ")
expmutm4_test <- subset(expmutm4_test, !duplicated(expmutm4_test[,1]))
summary(expmutm4_test)
#m5 S21
mutm5.mini <- as.data.frame(expmutm5merge[,1])
mutm5.mini <- cbind(mutm5.mini,total_reads=expmutm5merge$total_reads)
names(mutm5.mini)[names(mutm5.mini) == 'expmutm5merge[, 1]'] <- 'RECONSTRUCTED_SEQ'
mutm5.mini <- unique(mutm5.mini)
expmutm5_test<- merge(expmutm5,mutm5.mini, by="RECONSTRUCTED_SEQ")
expmutm5_test <- subset(expmutm5_test, !duplicated(expmutm5_test[,1]))
summary(expmutm5_test)
#m4 S22
wtm4.mini <- as.data.frame(expwtm4merge[,1])
wtm4.mini <- cbind(wtm4.mini,total_reads=expwtm4merge$total_reads)
names(wtm4.mini)[names(wtm4.mini) == 'expwtm4merge[, 1]'] <- 'RECONSTRUCTED_SEQ'
wtm4.mini <- unique(wtm4.mini)
expwtm4_test<- merge(expwtm4,wtm4.mini, by="RECONSTRUCTED_SEQ")
expwtm4_test <- subset(expwtm4_test, !duplicated(expwtm4_test[,1]))
summary(expwtm4_test)
#m5 S22
wtm5.mini <- as.data.frame(expwtm5merge[,1])
wtm5.mini <- cbind(wtm5.mini,total_reads=expwtm5merge$total_reads)
names(wtm5.mini)[names(wtm5.mini) == 'expwtm5merge[, 1]'] <- 'RECONSTRUCTED_SEQ'
wtm5.mini <- unique(wtm5.mini)
expwtm5_test<- merge(expwtm5,wtm5.mini, by="RECONSTRUCTED_SEQ")
expwtm5_test <- subset(expwtm5_test, !duplicated(expwtm5_test[,1]))
summary(expwtm5_test)
#m4 control (OregonR)
conm4.mini <- as.data.frame(conm4merge[,1])
conm4.mini <- cbind(conm4.mini,total_reads=conm4merge$total_reads)
names(conm4.mini)[names(conm4.mini) == 'conm4merge[, 1]'] <- 'RECONSTRUCTED_SEQ'
conm4.mini <- unique(conm4.mini)
conm4_test<- merge(conm4,conm4.mini, by="RECONSTRUCTED_SEQ")
conm4_test <- subset(conm4_test, !duplicated(conm4_test[,1]))
summary(conm4_test)
#m5 control (OregonR)
conm5.mini <- as.data.frame(conm5merge[,1])
conm5.mini <- cbind(conm5.mini,total_reads=conm5merge$total_reads)
names(conm5.mini)[names(conm5.mini) == 'conm5merge[, 1]'] <- 'RECONSTRUCTED_SEQ'
conm5.mini <- unique(conm5.mini)
conm5_test<- merge(conm5,conm5.mini, by="RECONSTRUCTED_SEQ")
conm5_test <- subset(conm5_test, !duplicated(conm5_test[,1]))
summary(conm5_test)

# remove any junctions with less thatn 30 reads
expmutm4_30 <- subset(expmutm4_test,total_reads>29)
summary(expmutm4_30)
expmutm5_30 <- subset(expmutm5_test,total_reads>29)
summary(expmutm5_30)
expwtm4_30 <- subset(expwtm4_test,total_reads>29)
summary(expwtm4_30)
expwtm5_30 <- subset(expwtm5_test,total_reads>29)
summary(expwtm5_30)
conm4_30 <- subset(conm4_test,total_reads>29)
summary(conm4_30)
conm5_30 <- subset(conm5_test,total_reads>29)
summary(conm5_30)

# Time to remove the control sequences from the experimentals
#M4 first
#m4 s21
toMatchm4 <- conm4_30$RECONSTRUCTED_SEQ
expmutm4_30$is.con<- regexpr(paste(toMatchm4,collapse="|"), expmutm4_30$RECONSTRUCTED_SEQ)
head(expmutm4_30)
con.expmutm4_30 <- subset(expmutm4_30, is.con=="-1")
#m4 s22
expwtm4_30$is.con<- regexpr(paste(toMatchm4,collapse="|"), expwtm4_30$RECONSTRUCTED_SEQ)
head(expwtm4_30)
con.expwtm4_30 <- subset(expwtm4_30, is.con=="-1")
#m5 s21
toMatchm5 <- conm5_30$RECONSTRUCTED_SEQ
expmutm5_30$is.con<- regexpr(paste(toMatchm5,collapse="|"), expmutm5_30$RECONSTRUCTED_SEQ)
head(expmutm5_30)
con.expmutm5_30 <- subset(expmutm5_30, is.con=="-1")
#m5 s22
expwtm5_30$is.con<- regexpr(paste(toMatchm5,collapse="|"), expwtm5_30$RECONSTRUCTED_SEQ)
head(expwtm5_30)
con.expwtm5_30 <- subset(expwtm5_30, is.con=="-1")

# combine all of them into one data table
con.master30 <- rbind(con.expmutm4_30,con.expmutm5_30,con.expwtm4_30,con.expwtm5_30)
write.csv(con.master30,"control_removed_data.csv")

m30 <- ggplot(con.master30, aes(x = -TOTAL_DELETION, colour = GENOTYPE))
m30 + geom_density(alpha=0.01, size=1) + geom_bar(stat='density', alpha=0.01,position="identity")+ facet_grid(PLASMID~ .)


#============================================================================================================================================
# August 1, 2017

#Determine percent of total and inaccurate reads for each junction
# To do this, we need to also include the exact match, which is removed in our con...._30 tables
# So we will add it back to that table for a little bit
summary(expmutm4_30)
exact.expmutm4_30 <- subset(expmutm4_30, CLASS=="exact")
exact.expmutm5_30 <- subset(expmutm5_30, CLASS=="exact")
exact.expwtm4_30 <- subset(expwtm4_30, CLASS=="exact")
exact.expwtm5_30 <- subset(expwtm5_30, CLASS=="exact")
exact.con.master30 <- rbind(con.master30,exact.expmutm4_30,exact.expmutm5_30,exact.expwtm4_30,exact.expwtm5_30)
write.csv(exact.con.master30,"control_removed_data_with_exact.csv")
exact.con.master30 <- read.csv("control_removed_data_with_exact.csv")
mut <- subset(exact.con.master30,GENOTYPE=="mutant")
wt <- subset(exact.con.master30,GENOTYPE=="wt")
mutm4 <- subset(mut, PLASMID=="M4")
mutm5 <- subset(mut,PLASMID=="M5")
wtm4 <- subset(wt,PLASMID=="M4")
wtm5 <- subset(wt,PLASMID=="M5")
mutm4$percent <- mutm4$X._OF_READS/sum(mutm4$X._OF_READS)*100
mutm5$percent <- mutm5$X._OF_READS/sum(mutm5$X._OF_READS)*100
wtm4$percent <- wtm4$X._OF_READS/sum(wtm4$X._OF_READS)*100
wtm5$percent <- wtm5$X._OF_READS/sum(wtm5$X._OF_READS)*100
sum(mutm4$percent)
sum(mutm5$percent)
sum(wtm4$percent)
sum(wtm5$percent)
# only calculate for percent of inaccurate reads
# m4 s21
mutm4.nexact <- subset(mutm4, CLASS!="exact")
mutm4.exact <- subset(mutm4,CLASS=="exact")
mutm4.exact$percent_inaccurate <- "NA"
mutm4.nexact$percent_inaccurate <- mutm4.nexact$total_reads/sum(mutm4.nexact$total_reads)*100
mutm4 <- rbind(mutm4.exact,mutm4.nexact)
# m5 s21
mutm5.nexact <- subset(mutm5, CLASS!="exact")
mutm5.exact <- subset(mutm5,CLASS=="exact")
mutm5.exact$percent_inaccurate <- "NA"
mutm5.nexact$percent_inaccurate <- mutm5.nexact$total_reads/sum(mutm5.nexact$total_reads)*100
mutm5 <- rbind(mutm5.exact,mutm5.nexact)
# m4 s22
wtm4.nexact <- subset(wtm4, CLASS!="exact")
wtm4.exact <- subset(wtm4,CLASS=="exact")
wtm4.exact$percent_inaccurate <- "NA"
wtm4.nexact$percent_inaccurate <- wtm4.nexact$total_reads/sum(wtm4.nexact$total_reads)*100
wtm4 <- rbind(wtm4.exact,wtm4.nexact)
# m5 s22
wtm5.nexact <- subset(wtm5, CLASS!="exact")
wtm5.exact <- subset(wtm5,CLASS=="exact")
wtm5.exact$percent_inaccurate <- "NA"
wtm5.nexact$percent_inaccurate <- wtm5.nexact$total_reads/sum(wtm5.nexact$total_reads)*100
wtm5 <- rbind(wtm5.exact,wtm5.nexact)
#combine them all back together again and save new data table
p.master30 <- rbind(mutm4,mutm5,wtm4,wtm5)
write.csv(p.master30,"control_removed_data_percent.csv")

# take a look at different types of classes of junctions
tbls2 <- aggregate(as.numeric(percent_inaccurate)~CLASS+PLASMID+GENOTYPE,data=p.master30, sum)
names(tbls2)[names(tbls2) == 'as.numeric(percent_inaccurate)'] <- 'percent_inaccurate'

f1 <- ggplot(tbls2, aes(fill=CLASS))+ 
  geom_bar(aes(x=PLASMID, y=percent_inaccurate), position="stack", stat="identity", colour="black") +
  facet_wrap(~GENOTYPE)+
  theme_classic(base_size = 9)+
  guides(fill=FALSE)+
  scale_fill_grey(start=0.0, end=1)+
  scale_y_continuous(name = "% Inaccurate Reads",
                     labels=c("0","25","50","75","100"),
                     expand=c(0,0),
                     limits=c(0,101),
                     breaks=c(0,25,50,75,100))+
  scale_x_discrete(name = "Plasmid")+
  theme(axis.text.x=element_text(size=9, face="bold"),
        axis.text.y=element_text(size=9, face="bold"),
        axis.title.x=element_text(size=10, face="bold"),
        axis.title.y=element_text(size=10, face="bold"))
f1
# not very useful, because we need to actually break down complex into either insertion or deletion, and break down deletions into MHJ or ABJ


#-------------Get deletion sequences for Python Program
# some alignments characterized as "complex" events are just misaligned deletion events. I need to get these out and run them through the python program
exp.30 <- subset(p.master30, SAMPLE=="exp") # repeat from above line just in case
complex <- subset(exp.30, CLASS=="complex")
# i think think I can begin to seperate complex indels from deletions by removing any complex sequences with a net possitive insertion. those would be true insertions
complex$net <- complex$TOTAL_DELETION+complex$INSERTION_LENGTH
comp.true <- subset(complex, net<1)
comp.false <- subset(complex,net>0)
comp.false$CLASS <- "insertion"
mutcomp <- subset(comp.true, GENOTYPE=="mutant")
wtcomp <- subset(comp.true, GENOTYPE=="wt")
mutm4comp <- subset(mutcomp, PLASMID=="M4")
mutm5comp <- subset(mutcomp, PLASMID=="M5")
wtm4comp <- subset(wtcomp, PLASMID=="M4")
wtm5comp <- subset(wtcomp, PLASMID=="M5")

write.csv(mutm4comp, "mut_m4_complex_del.csv")
write.csv(mutm5comp, "mut_m5_complex_del.csv")
write.csv(wtm4comp, "wt_m4_complex_del.csv")
write.csv(wtm5comp, "wt_m5_complex_del.csv")

# I went through "iw7_complex_del.csv" spreadsheet and manually checked if deletion or not and created a collumn "is.complex" and saved. If complex wrote "-C"
# time to upload, seperate based on "is.complex" and then run through SD-MMEJ program
mutm4curate <- read.csv("mut_m4_complex_del-curate2.csv")
mutm5curate <- read.csv("mut_m5_complex_del-curate2.csv")
wtm4curate <- read.csv("wt_m4_complex_del-curate2.csv")
wtm5curate <- read.csv("wt_m5_complex_del-curate2.csv")

#exclude complex
idc <- subset(mutm4curate, is.complex!="-C") 
m1dc <- subset(mutm5curate, is.complex!="-C")
m2dc <- subset(wtm4curate, is.complex!="-C")
m3dc <- subset(wtm5curate, is.complex!="-C")
#reclasify as deletion
idc$CLASS <- "deletion"
m1dc$CLASS <- "deletion"
m2dc$CLASS <- "deletion"
m3dc$CLASS <- "deletion"
#add to deletion subset
del <- subset(p.master30, CLASS=="deletion")
del2 <- rbind(del[,2:30],idc[,3:31],m1dc[,3:31],m2dc[,3:31],m3dc[,3:31])

# we need to also reclassify the complex ones as insertions and add to the insertion subset
iic <- subset(mutm4curate, is.complex=="-C") 
m1ic <- subset(mutm5curate, is.complex=="-C")
m2ic <- subset(wtm4curate, is.complex=="-C")
m3ic <- subset(wtm5curate, is.complex=="-C")
#reclasify as deletion
iic$CLASS <- "insertion"
m1ic$CLASS <- "insertion"
m2ic$CLASS <- "insertion"
m3ic$CLASS <- "insertion"
ins <- subset(p.master30, CLASS=="insertion")
ins2 <- rbind(ins[,2:30],iic[,3:31],m1ic[,3:31],m2ic[,3:31],m3ic[,3:31])
fixed <- rbind(ins2,del2,comp.false[,2:30])
sum(as.numeric(fixed$percent_inaccurate))
exact <- subset(p.master30, CLASS=="exact")
fixed_exact <- rbind(fixed,exact[,2:30])
sum(as.numeric(fixed_exact$percent))
write.csv(fixed_exact,"combined_curated.csv")
# this doesnt take into account the ones

#============================================================================================================================
# August 2, 2017

# Get the deletion sequences for python
del <- subset(fixed_exact, CLASS=="deletion")
mutdel <- subset(del, GENOTYPE=="mutant")
wtdel <- subset(del, GENOTYPE=="wt")
mutm4del <- subset(mutdel, PLASMID=="M4")
mutm5del <- subset(mutdel, PLASMID=="M5")
wtm4del <- subset(wtdel, PLASMID=="M4")
wtm5del <- subset(wtdel, PLASMID=="M5")

mutm4del.seq <- as.data.frame(mutm4del$ALIGNED_SEQ)
mutm5del.seq <- as.data.frame(mutm5del$ALIGNED_SEQ)
wtm4del.seq <- as.data.frame(wtm4del$ALIGNED_SEQ)
wtm5del.seq <- as.data.frame(wtm5del$ALIGNED_SEQ)
#save the sequences to put through Amy's SD-MMEJ program
write.table(mutm4del.seq, "mut_m4_del_for_python.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(mutm5del.seq, "mut_m5_del_for_python.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(wtm4del.seq, "wt_m4_del_for_python.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(wtm5del.seq, "wt_m5_del_for_python.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

# Run sequences through program, and combine into excel sheet
# Combined deletions from the SD-MMEJ program output into a file with the sequence, deletion length, left/right deletion index, and consistency
p <- read.csv("m4_m5_amplicon_python_consistency.csv")
p$SEQ <- toupper(p$SEQ)
# seperate based on plasmid and genotype
mutp <- subset(p,GENOTYPE=="mutant")
wtp <- subset(p,GENOTYPE=="wt")
mutm4p <- subset(mutp, PLASMID=="M4")
mutm5p <- subset(mutp, PLASMID=="M5")
wtm4p <- subset(wtp, PLASMID=="M4")
wtm5p <- subset(wtp, PLASMID=="M5")

# Need to merge the consistnecy data with the data from Steves analysis to give us number of reads etc.

#merge sequences from SD-MMEJ python program with the above del datasets
# Now try the above consistency analysis with the Percent inaccurate events
mutm4p.merge <- merge(mutm4del, mutm4p, by.x="RECONSTRUCTED_SEQ", by.y="SEQ", all=FALSE)
mutm5p.merge <- merge(mutm5del, mutm5p, by.x="RECONSTRUCTED_SEQ", by.y="SEQ", all=FALSE)
wtm4p.merge <- merge(wtm4del, wtm4p, by.x="RECONSTRUCTED_SEQ", by.y="SEQ", all=FALSE)
wtm5p.merge <- merge(wtm5del, wtm5p, by.x="RECONSTRUCTED_SEQ", by.y="SEQ", all=FALSE)

#save these. i forget when...but you will want to come back to these
write.csv(mutm4p.merge,"mut_m4_del_merge_table.csv")
write.csv(mutm5p.merge,"mut_m5_del_merge_table.csv")
write.csv(wtm4p.merge,"wt_m4_del_merge_table.csv")
write.csv(wtm5p.merge,"wt_m5_del_merge_table.csv")

# now that we have read information merged with SD-mmej consistency we can start to figure out % SD-MMEJ consistent
mutm4p.merge <- read.csv("mut_m4_del_merge_table.csv")
mutm5p.merge <- read.csv("mut_m5_del_merge_table.csv")
wtm4p.merge <- read.csv("wt_m4_del_merge_table.csv")
wtm5p.merge <- read.csv("wt_m5_del_merge_table.csv")
# Need to 


#------------ M4 S21
# YAY it works! Now figure out how to do weighted counts or something
max.ip <- sum(as.numeric(mutm4p.merge$percent_inaccurate))
max.ip
# total percentage of innacurate reads are deletion jucntions = [1] 76.91989
abj.i <- subset(mutm4p.merge, REPAIR_TYPE=="ABJ")
abj.ip <- sum(as.numeric(abj.i$percent_inaccurate))
abj.ip
# percent of inaccurate reads that are ABJ = 22.66611
rel.abj <- abj.ip/max.ip*100 
rel.abj
# 29.46717% of Deletion Junctions are ABJ
mhj.i <- subset(mutm4p.merge, REPAIR_TYPE=="MHJ")
mhj.ip <- sum(as.numeric(mhj.i$percent_inaccurate))
mhj.ip
# percent of inaccurate reads that are MHJ = 54.25378
rel.mhj <- mhj.ip/max.ip*100 
rel.mhj
# 70.53283% of Deletion junctions are MHJ
abj.c <- sum(as.numeric(abj.i[which(abj.i[,37]=="TRUE"),30]))
abj.c # 20.55475
abj.c/abj.ip*100
#[1]  90.68493 Of ABJ are consistent
mhj.c <- sum(as.numeric(mhj.i[which(mhj.i[,37]=="TRUE"),30]))
mhj.c # 30.81143
mhj.c/mhj.ip*100
#[1] 56.7913% of MHJ are consistent
#total consistent deletions
tot.c <- sum(as.numeric(mutm4p.merge[which(mutm4p.merge[,37]=="TRUE"),30]))
tot.c # 51.36618% of inaccurate Deletions reads that are consistent
tot.c/max.ip*100
# [1] 66.77879% of Deletions are consistent

#------------ M5 S21
# YAY it works! Now figure out how to do weighted counts or something
max.ip <- sum(as.numeric(mutm5p.merge$percent_inaccurate))
max.ip
# total percentage of innacurate reads are deletion jucntions = [1] 66.88356
abj.i <- subset(mutm5p.merge, REPAIR_TYPE=="ABJ")
abj.ip <- sum(as.numeric(abj.i$percent_inaccurate))
abj.ip
# percent of inaccurate reads that are ABJ = 23.5274
rel.abj <- abj.ip/max.ip*100 
rel.abj
# 35.17665% of Deletion Junctions are ABJ
mhj.i <- subset(mutm5p.merge, REPAIR_TYPE=="MHJ")
mhj.ip <- sum(as.numeric(mhj.i$percent_inaccurate))
mhj.ip
# percent of inaccurate reads that are MHJ = 43.35616
rel.mhj <- mhj.ip/max.ip*100 
rel.mhj
# 64.82335% of Deletion junctions are MHJ
abj.c <- sum(as.numeric(abj.i[which(abj.i[,37]=="TRUE"),30]))
abj.c # 22.24886
abj.c/abj.ip*100
#[1]  94.56574 Of ABJ are consistent
mhj.c <- sum(as.numeric(mhj.i[which(mhj.i[,37]=="TRUE"),30]))
mhj.c # 31.74658
mhj.c/mhj.ip*100
#[1] 73.22275% of MHJ are consistent
#total consistent deletions
tot.c <- sum(as.numeric(mutm5p.merge[which(mutm5p.merge[,37]=="TRUE"),30]))
tot.c # 53.99543% of inaccurate Deletions reads that are consistent
tot.c/max.ip*100
# [1] 80.7305% of Deletions are consistent

#------------ M4 S22 (WT)
# YAY it works! Now figure out how to do weighted counts or something
max.ip <- sum(as.numeric(wtm4p.merge$percent_inaccurate))
max.ip
# total percentage of innacurate reads are deletion jucntions = [1] 58.45993
abj.i <- subset(wtm4p.merge, REPAIR_TYPE=="ABJ")
abj.ip <- sum(as.numeric(abj.i$total_reads))
abj.ip
# percent of inaccurate reads that are ABJ = 19.03485
rel.abj <- abj.ip/max.ip*100 
rel.abj
# 31.77558% of Deletion Junctions are ABJ
mhj.i <- subset(wtm4p.merge, REPAIR_TYPE=="MHJ")
mhj.ip <- sum(as.numeric(mhj.i$total_reads))
mhj.ip
# percent of inaccurate reads that are MHJ = 39.42508
rel.mhj <- mhj.ip/max.ip*100 
rel.mhj
# 68.22442% of Deletion junctions are MHJ
abj.c <- sum(as.numeric(abj.i[which(abj.i[,37]=="TRUE"),27]))
abj.c # 17.30712
abj.c/abj.ip*100
#[1]  90.09044 Of ABJ are consistent
mhj.c <- sum(as.numeric(mhj.i[which(mhj.i[,37]=="TRUE"),27]))
mhj.c # 27.50968
mhj.c/mhj.ip*100
#[1] 64.67033% of MHJ are consistent
#total consistent deletions
tot.c <- sum(as.numeric(wtm4p.merge[which(wtm4p.merge[,37]=="TRUE"),27]))
tot.c # 44.8168% of inaccurate Deletions reads that are consistent
tot.c/max.ip*100
# [1] 72.74772% of Deletions are consistent

#------------ M5 S22 (WT)
# YAY it works! Now figure out how to do weighted counts or something
max.ip <- sum(as.numeric(wtm5p.merge$total_reads))
max.ip
# total percentage of innacurate reads are deletion jucntions = [1] 67.82658
abj.i <- subset(wtm5p.merge, REPAIR_TYPE=="ABJ")
abj.ip <- sum(as.numeric(abj.i$total_reads))
abj.ip
# percent of inaccurate reads that are ABJ = 29.94866
rel.abj <- abj.ip/max.ip*100 
rel.abj
# 34.74178% of Deletion Junctions are ABJ
mhj.i <- subset(wtm5p.merge, REPAIR_TYPE=="MHJ")
mhj.ip <- sum(as.numeric(mhj.i$total_reads))
mhj.ip
# percent of inaccurate reads that are MHJ = 37.87792
rel.mhj <- mhj.ip/max.ip*100 
rel.mhj
# 65.25822% of Deletion junctions are MHJ
abj.c <- sum(as.numeric(abj.i[which(abj.i[,37]=="TRUE"),27]))
abj.c # 28.3514
abj.c/abj.ip*100
#[1]  93.45114 Of ABJ are consistent
mhj.c <- sum(as.numeric(mhj.i[which(mhj.i[,37]=="TRUE"),27]))
mhj.c # 34.62635
mhj.c/mhj.ip*100
#[1] 72.88323% of MHJ are consistent
#total consistent deletions
tot.c <- sum(as.numeric(wtm5p.merge[which(wtm5p.merge[,37]=="TRUE"),27]))
tot.c # 62.97775% of inaccurate Deletions reads that are consistent
tot.c/max.ip*100
# [1] 80.02889% of Deletions are consistent


#=================================================================================================================
# August 3, 2017
# analyze insertion junctions

#------Insertion Analysis
#Need to get the insertion sequences to run through my R program "INSERTION_PRORGAM_20160725"
# in theory you could just tell the program to pick out this collumn....but i swear to god if you try and do that and it fucks up, dont ask me for help.
# get RECONSTRUCTED_SEQ to run through python SD-MMEJ program
ins <- subset(fixed_exact, CLASS=="insertion")
mut.i <- subset(ins, GENOTYPE=="mutant")
wt.i <- subset(ins, GENOTYPE=="wt")
mutm4.i <- subset(mut.i, PLASMID=="M4")
mutm5.i <- subset(mut.i, PLASMID=="M5")
wtm4.i <- subset(wt.i, PLASMID=="M4")
wtm5.i <- subset(wt.i, PLASMID=="M5")

mutm4.i.seq <- as.data.frame(mutm4.i$RECONSTRUCTED_SEQ)
names(mutm4.i.seq) <- "RECONSTRUCTED_SEQ"
mutm5.i.seq <- as.data.frame(mutm5.i$RECONSTRUCTED_SEQ)
names(mutm5.i.seq) <- "RECONSTRUCTED_SEQ"
wtm4.i.seq <- as.data.frame(wtm4.i$RECONSTRUCTED_SEQ)
names(wtm4.i.seq) <- "RECONSTRUCTED_SEQ"
wtm5.i.seq <- as.data.frame(wtm5.i$RECONSTRUCTED_SEQ)
names(wtm5.i.seq) <- "RECONSTRUCTED_SEQ"

# save to run through R insertion analysis program. This is a seperate script designed to only need to change certain parameters at the top and then run through. Make sure you have an version that you never go to unless you fuck up something. I dont know how i managed to make this script in the first place, and i dont know how i will be able to fix it.
write.csv(mutm4.i.seq, "mut_m4_ins_for_R_script.csv",row.names=FALSE)
write.csv(mutm5.i.seq, "mut_m5_ins_for_R_script.csv",row.names=FALSE)
write.csv(wtm4.i.seq, "wt_m4_ins_for_R_script.csv",row.names=FALSE)
write.csv(wtm5.i.seq, "wt_m5_ins_for_R_script.csv",row.names=FALSE)


#-----------------------------------------------------------------------------
# Ran the through the insertion program for each plasmid. I had to go through manually to change some indexes for the start and end of repeat motifs. Also to change some consistency. Then i reran through the end of the insertion program to remake the tables. They are labeled with a "2" at the end. This step are in the individual R scripts saved in the insertion_analysis folder

# I made an R script that does analysis of insertions that are single step insertions, and determine if they are SD-MMEJ consistent
# Now add these data to the deletion data and do some fucking stats!
setwd("R:/Varandt/SD-MMEJ/07132017_amplicon_sequence")
m <- read.csv("combined_curated.csv")
a <- read.csv("mut_m4_insertion_consistency3.csv")
b <- read.csv("mut_m5_insertion_consistency4.csv")
c <- read.csv("wt_m4_insertion_consistency4.csv")
d <- read.csv("wt_m5_insertion_consistency3.csv")

# The deletion boundaries are totally off on the insertions, I need to take into account all of the insertions and shit.
a$RIGHT_DEL <- 165-(nchar(as.character(a$RECONSTRUCTED_SEQ))-a$right_del+1)
b$RIGHT_DEL <- 165-(nchar(as.character(b$RECONSTRUCTED_SEQ))-b$right_del+1)
c$RIGHT_DEL <- 165-(nchar(as.character(c$RECONSTRUCTED_SEQ))-c$right_del+1)
d$RIGHT_DEL <- 165-(nchar(as.character(d$RECONSTRUCTED_SEQ))-d$right_del+1)

ins <- subset(m, CLASS=="insertion")
mut.i <- subset(ins, GENOTYPE=="mutant")
wt.i <- subset(ins, GENOTYPE=="wt")
mutm4.i <- subset(mut.i,PLASMID=="M4")
mutm5.i <- subset(mut.i,PLASMID=="M5")
wtm4.i <- subset(wt.i,PLASMID=="M4")
wtm5.i <- subset(wt.i,PLASMID=="M5")


a$plasmid <- "M4"
b$plasmid <- "M5"
c$plasmid <- "M4"
d$plasmid <- "M5"

# Need to remove the insertion sequences where it says TRUE and FALSE...they just need to all say TRUE, remove the row that says FALSE
asort <-a[order(a$consistency, decreasing=TRUE),]
ashort <- asort[!duplicated(asort[,8]),]
ashort2 <- ashort[,7:9]
ashort2$inserted_seq <- ashort$insertion
ashort2$right_del <- ashort[,18]
mutm4a <- merge(mutm4.i, ashort2, by="RECONSTRUCTED_SEQ")
# Calculate percent consistency
max.iip <- sum(mutm4a$percent_inaccurate)
max.iip
# total percentage of innacurate reads are insertions or indels = 23.08011%
ins.c <- sum(mutm4a[which(mutm4a[,31]=="TRUE"),30])
ins.c # 8.186711% of insertion reads are consistent
ins.c/max.iip*100
#[1] 35.47085
# Percent of one step SD-MMEJ consistent insertions is 35.47085%

# M5 S21
# Need to remove the insertion sequences where it says TRUE and FALSE...they just need to all say TRUE, remove the row that says FALSE
bsort <-b[order(b$consistency, decreasing=TRUE),]
bshort <- bsort[!duplicated(bsort[,8]),]
bshort2 <- bshort[,7:9]
bshort2$inserted_seq <- bshort$insertion
bshort2$right_del <- bshort[,18]
mutm5a <- merge(mutm5.i, bshort2, by="RECONSTRUCTED_SEQ")
# Calculate percent consistency
max.m1ip <- sum(mutm5a$percent_inaccurate)
max.m1ip
# total percentage of innacurate reads are insertions or indels = 33.11644%
ins.c <- sum(mutm5a[which(mutm5a[,31]=="TRUE"),30])
ins.c # 15.1484
ins.c/max.m1ip*100
#[1] 45.74285
# Percent of one step SD-MMEJ consistent insertions is 45.74285%

# M4 S22
# Need to remove the insertion sequences where it says TRUE and FALSE...they just need to all say TRUE, remove the row that says FALSE
csort <-c[order(c$consistency, decreasing=TRUE),]
cshort <- csort[!duplicated(csort[,8]),]
cshort2 <- cshort[,7:9]
cshort2$inserted_seq <- cshort$insertion
cshort2$right_del <- cshort[,18]
wtm4a <- merge(wtm4.i, cshort2, by="RECONSTRUCTED_SEQ")
# Calculate percent consistency
max.m2ip <- sum(wtm4a$percent_inaccurate)
max.m2ip
# total percentage of innacurate reads are insertions or indels = 35.01421%
ins.c <- sum(wtm4a[which(wtm4a[,31]=="TRUE"),30])
ins.c # 16.02923
ins.c/max.m2ip*100
#[1] 45.77922
# Percent of one step SD-MMEJ consistent insertions is 45.77922%

# M5 S22
# Need to remove the insertion sequences where it says TRUE and FALSE...they just need to all say TRUE, remove the row that says FALSE
dsort <-d[order(d$consistency, decreasing=TRUE),]
dshort <- dsort[!duplicated(dsort[,8]),]
dshort2 <- dshort[,7:9]
dshort2$inserted_seq <- dshort$insertion
dshort2$right_del <- dshort[,18]
wtm5a <- merge(wtm5.i, dshort2, by="RECONSTRUCTED_SEQ")
# Calculate percent consistency
max.m3ip <- sum(wtm5a$percent_inaccurate)
max.m3ip
# total percentage of innacurate reads are insertions or indels = 23.27515%
ins.c <- sum(wtm5a[which(wtm5a[,31]=="TRUE"),30])
ins.c # 11.13882
ins.c/max.m3ip*100
#[1] 47.85714
# Percent of one step SD-MMEJ consistent insertions is 47.85714%


#===============Need to combine deletion and insertion consistency tables to create 1 table that has all of the consistency on it with accurate deletion boundaries.

del <- read.csv("m4_m5_del_merge_table_consistency.csv") # this table was made by combining "mut_m4_del_merge_table.csv", "mut_m5_del_merge_table.csv", "wt_m4_del_merge_table.csv", and "wt_m5_del_merge_table.csv" in Excel and removing the row numbers column, along with the second set of PLASMID and GENOTYPE columns. they were unnecessary 

#Need to combine the del and insertion tables into one that includes consistency
ins <- rbind(mutm4a,mutm5a,wtm4a,wtm5a)
ins1 <- ins[,3:21]
ins1$inserted_seq <- ins$inserted_seq
ins1$RECONSTRUCTED_SEQ <- ins$RECONSTRUCTED_SEQ
ins1$X._OF_READS <- ins$X._OF_READS
ins1$total_reads <- ins$total_reads
ins1$MICROHOMOLOGY <- ins$MICROHOMOLOGY 
ins1$MH_Length <- ins$MH_Length
ins1$NUMBER_OF_ALIGNMENTS  <- ins$NUMBER_OF_ALIGNMENTS
ins1$ALIGNED_SEQ <- ins$ALIGNED_SEQ
ins1$percent <- ins$percent
ins1$percent_inaccurate <- ins$percent_inaccurate
ins1$REPAIR_TYPE <- "InDel"
ins1$left_del <- ins$left_del 
ins1$right_del <- ins$right_del
ins1$CONSISTENCY <- ins$consistency

del1 <- cbind(del[,2:20])
del1$inserted_seq <- NA
del1$RECONSTRUCTED_SEQ <- del$RECONSTRUCTED_SEQ
del1$X._OF_READS <- del$X._OF_READS
del1$total_reads <- del$total_reads
del1$MICROHOMOLOGY <- del$MICROHOMOLOGY 
del1$MH_Length <- del$MH_Length
del1$NUMBER_OF_ALIGNMENTS  <- del$NUMBER_OF_ALIGNMENTS
del1$ALIGNED_SEQ <- del$ALIGNED_SEQ
del1$percent <- del$percent
del1$percent_inaccurate <- del$percent_inaccurate
del1$REPAIR_TYPE <- del$REPAIR_TYPE
del1 <- cbind(del1, del[,32:34])
names(del1) <- names(ins1)
all <- rbind(ins1,del1)
sum(all$percent_inaccurate)
write.csv(all,"m4_m5_amplicon_SD-MMEJ_consistency2.csv")

#seperate based on plasmid
mut <- subset(all, GENOTYPE=="mutant")
wt <- subset(all, GENOTYPE=="wt")
mutm4 <- subset(mut,PLASMID=="M4")
mutm5 <- subset(mut,PLASMID=="M5")
wtm4 <- subset(wt,PLASMID=="M4")
wtm5 <- subset(wt,PLASMID=="M5")

#calculate percent consistent for ALL junction types
mutm4.c <- sum(mutm4[which(mutm4[,33]=="TRUE"),29])
mutm4.c # 59.55289
mutm5.c <- sum(mutm5[which(mutm5[,33]=="TRUE"),29])
mutm5.c # 69.14384
wtm4.c <- sum(wtm4[which(wtm4[,33]=="TRUE"),29])
wtm4.c # 63.30491
wtm5.c <- sum(wtm5[which(wtm5[,33]=="TRUE"),29])
wtm5.c # 72.54087





#------------------------------% of Sanger sequences in Amplicon-----------------------
##I should probably start to compare these to the sanger data!
#setwd("R:/Varandt/SD-MMEJ/20160616_amplicon_sequence/Fixed_MHanalysis_11_20_16")
#amp <- read.csv("amplicon_SD-MMEJ_consistency.csv")
#aiw7 <- subset(amp, PLASMID=="Iw7")
#am1 <- subset(amp, PLASMID=="M1")
#am2 <- subset(amp, PLASMID=="M2")
#am3 <- subset(amp, PLASMID=="M3")
#sd<-read.csv("sanger_del_bound.csv")
#summary(sd)
##need to remove all non alphabet characters from the aligned seq and make uppercase
#clean <- gsub('[0-9]+', '', sd$ALIGNED_SEQ)
#clean <- str_replace_all(clean, "[[:punct:]]", "")
#clean <- toupper(clean)
#sd$ALIGNED_SEQ2 <- as.factor(clean)
#siw7 <- subset(sd, PLASMID=="iw7")
#sm1 <- subset(sd, PLASMID=="m1")
#sm2 <- subset(sd, PLASMID=="m2")
#sm3 <- subset(sd, PLASMID=="m3")

## Need to determine the % of sanger jxns in our amplicon data
## remove duplicate sanger junctions
#single.siw7 <- subset(siw7,!duplicated(siw7$ALIGNED_SEQ2))
#toMatch <- single.siw7$ALIGNED_SEQ
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      aiw7$RECONSTRUCTED_SEQ, value=TRUE)))
##21/43 matches --> That is pretty low.
## lets see if we only look at Deletions or insertions
## Deletions first
#del.siw7 <- subset(single.siw7, REPAIR_TYPE != "InDel")
#toMatch <- del.siw7$ALIGNED_SEQ
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      aiw7$RECONSTRUCTED_SEQ, value=TRUE)))
##18/21 sanger deletions are found in amplicons
## Insertions
#indel.siw7 <- subset(single.siw7, REPAIR_TYPE== "InDel")
#toMatch <- indel.siw7$ALIGNED_SEQ2
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      aiw7$RECONSTRUCTED_SEQ, value=TRUE)))
## 2/22 insertions
## M1
## Need to determine the % of sanger jxns in our amplicon data
## remove duplicate sanger junctions
#single.sm1 <- subset(sm1,!duplicated(sm1$ALIGNED_SEQ2))
#toMatch <- single.sm1$ALIGNED_SEQ2
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      am1$RECONSTRUCTED_SEQ, value=TRUE)))
##21/39 matches --> That is pretty low.
## lets see if we only look at Deletions or insertions
## Deletions first
#del.sm1 <- subset(single.sm1, REPAIR_TYPE != "InDel")
#toMatch <- del.sm1$ALIGNED_SEQ
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      am1$RECONSTRUCTED_SEQ, value=TRUE)))
##13/14 sanger deletions are found in amplicons
## Insertions
#indel.sm1 <- subset(single.sm1, REPAIR_TYPE== "InDel")
#toMatch <- indel.sm1$ALIGNED_SEQ2
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      am1$RECONSTRUCTED_SEQ, value=TRUE)))
##8/25 insertions

## M2
## Need to determine the % of sanger jxns in our amplicon data
## remove duplicate sanger junctions
#single.sm2 <- subset(sm2,!duplicated(sm2$ALIGNED_SEQ2))
#toMatch <- single.sm2$ALIGNED_SEQ2
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      am2$RECONSTRUCTED_SEQ, value=TRUE)))
##27/47 matches --> That is pretty low.
## lets see if we only look at Deletions or insertions
## Deletions first
#del.sm2 <- subset(single.sm2, REPAIR_TYPE != "InDel")
#toMatch <- del.sm2$ALIGNED_SEQ
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      am2$RECONSTRUCTED_SEQ, value=TRUE)))
##18/19 sanger deletions are found in amplicons
## Insertions
#indel.sm2 <- subset(single.sm2, REPAIR_TYPE== "InDel")
#toMatch <- indel.sm2$ALIGNED_SEQ2
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      am2$RECONSTRUCTED_SEQ, value=TRUE)))
##9/28 insertions

## M3
## Need to determine the % of sanger jxns in our amplicon data
## remove duplicate sanger junctions
#single.sm3 <- subset(sm3,!duplicated(sm3$ALIGNED_SEQ2))
#toMatch <- single.sm3$ALIGNED_SEQ2
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      am3$RECONSTRUCTED_SEQ, value=TRUE)))
##16/40 matches --> That is pretty low.
## lets see if we only look at Deletions or insertions
## Deletions first
#del.sm3 <- subset(single.sm3, REPAIR_TYPE != "InDel")
#toMatch <- del.sm3$ALIGNED_SEQ
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      am3$RECONSTRUCTED_SEQ, value=TRUE)))
##10/10 sanger deletions are found in amplicons
## Insertions
#indel.sm3 <- subset(single.sm3, REPAIR_TYPE== "InDel")
#toMatch <- indel.sm3$ALIGNED_SEQ2
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      am3$RECONSTRUCTED_SEQ, value=TRUE)))
##6/30 insertions

## Need to try and see how many sanger sequences are in amplicon without the 30 read cut off
##I should probably start to compare these to the sanger data!
#setwd("R:/Varandt/SD-MMEJ/20160616_amplicon_sequence/Fixed_MHanalysis_11_20_16")
#amp <- read.csv("combined_fixed.csv")
#amp <- subset(amp, SAMPLE=="exp")
#aiw7 <- subset(amp, PLASMID=="Iw7")
#am1 <- subset(amp, PLASMID=="M1")
#am2 <- subset(amp, PLASMID=="M2")
#am3 <- subset(amp, PLASMID=="M3")
#sd<-read.csv("sanger_del_bound.csv")
#summary(sd)
##need to remove all non alphabet characters from the aligned seq and make uppercase
#clean <- gsub('[0-9]+', '', sd$ALIGNED_SEQ)
#clean <- str_replace_all(clean, "[[:punct:]]", "")
#clean <- toupper(clean)
#sd$ALIGNED_SEQ2 <- as.factor(clean)
#siw7 <- subset(sd, PLASMID=="iw7")
#sm1 <- subset(sd, PLASMID=="m1")
#sm2 <- subset(sd, PLASMID=="m2")
#sm3 <- subset(sd, PLASMID=="m3")

## Need to determine the % of sanger jxns in our amplicon data
## remove duplicate sanger junctions
#single.siw7 <- subset(siw7,!duplicated(siw7$ALIGNED_SEQ2))
#toMatch <- single.siw7$ALIGNED_SEQ
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      aiw7$RECONSTRUCTED_SEQ, value=TRUE)))
##31/43 matches --> That is pretty low.
## lets see if we only look at Deletions or insertions
## Deletions first
#del.siw7 <- subset(single.siw7, REPAIR_TYPE != "InDel")
#toMatch <- del.siw7$ALIGNED_SEQ
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      aiw7$RECONSTRUCTED_SEQ, value=TRUE)))
##18/21 sanger deletions are found in amplicons
## Insertions
#indel.siw7 <- subset(single.siw7, REPAIR_TYPE== "InDel")
#toMatch <- indel.siw7$ALIGNED_SEQ2
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      aiw7$RECONSTRUCTED_SEQ, value=TRUE)))
## 3/22 insertions
## M1
## Need to determine the % of sanger jxns in our amplicon data
## remove duplicate sanger junctions
#single.sm1 <- subset(sm1,!duplicated(sm1$ALIGNED_SEQ2))
#toMatch <- single.sm1$ALIGNED_SEQ2
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      am1$RECONSTRUCTED_SEQ, value=TRUE)))
##25/39 matches --> That is pretty low.
## lets see if we only look at Deletions or insertions
## Deletions first
#del.sm1 <- subset(single.sm1, REPAIR_TYPE != "InDel")
#toMatch <- del.sm1$ALIGNED_SEQ
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      am1$RECONSTRUCTED_SEQ, value=TRUE)))
##14/14 sanger deletions are found in amplicons
## Insertions
#indel.sm1 <- subset(single.sm1, REPAIR_TYPE== "InDel")
#toMatch <- indel.sm1$ALIGNED_SEQ2
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      am1$RECONSTRUCTED_SEQ, value=TRUE)))
##11/25 insertions

## M2
## Need to determine the % of sanger jxns in our amplicon data
## remove duplicate sanger junctions
#single.sm2 <- subset(sm2,!duplicated(sm2$ALIGNED_SEQ2))
#toMatch <- single.sm2$ALIGNED_SEQ2
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      am2$RECONSTRUCTED_SEQ, value=TRUE)))
##35/47 matches --> That is pretty low.
## lets see if we only look at Deletions or insertions
## Deletions first
#del.sm2 <- subset(single.sm2, REPAIR_TYPE != "InDel")
#toMatch <- del.sm2$ALIGNED_SEQ
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      am2$RECONSTRUCTED_SEQ, value=TRUE)))
##18/19 sanger deletions are found in amplicons
## Insertions
#indel.sm2 <- subset(single.sm2, REPAIR_TYPE== "InDel")
#toMatch <- indel.sm2$ALIGNED_SEQ2
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      am2$RECONSTRUCTED_SEQ, value=TRUE)))
##13/28 insertions

## M3
## Need to determine the % of sanger jxns in our amplicon data
## remove duplicate sanger junctions
#single.sm3 <- subset(sm3,!duplicated(sm3$ALIGNED_SEQ2))
#toMatch <- single.sm3$ALIGNED_SEQ2
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      am3$RECONSTRUCTED_SEQ, value=TRUE)))
##24/40 matches --> That is pretty low.
## lets see if we only look at Deletions or insertions
## Deletions first
#del.sm3 <- subset(single.sm3, REPAIR_TYPE != "InDel")
#toMatch <- del.sm3$ALIGNED_SEQ
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      am3$RECONSTRUCTED_SEQ, value=TRUE)))
##10/10 sanger deletions are found in amplicons
## Insertions
#indel.sm3 <- subset(single.sm3, REPAIR_TYPE== "InDel")
#toMatch <- indel.sm3$ALIGNED_SEQ2
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      am3$RECONSTRUCTED_SEQ, value=TRUE)))
##13/30 insertions


#single.sd <- subset(sd,!duplicated(sd$SEQ))
#head(single.sd)
##create "seq" without "---"
## to do this, i need to remove the '----' from the aligned sequences that were placed into python
#single.sd$seq <- gsub("-", "",single.sd$SEQ)
#siw7 <- subset(sd, PLASMID=="iw7")
#siw7 <- subset(single.sd, PLASMID=="Iw7")
##also need to make sure i am only looking at deletions for the amplicon
#iw7.del <- subset(s.ap.thirty, PLASMID=="Iw7")
#toMatch <- siw7$seq
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      iw7.del$RECONSTRUCTED_SEQ, value=TRUE)))
##14/21
#matches.all <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                          s.iw7$RECONSTRUCTED_SEQ, value=TRUE)))
##25/21? --> but this is for ALL jxns, without a cuttoff. Time to apply the cuttoff then redo it
#s.iw7.30 <- subset(s.iw7, X._OF_READS>29)
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      s.iw7.30$RECONSTRUCTED_SEQ, value=TRUE)))
##20/21! MUCH BETTER! This means that there are several "complex" jxns that are just deletions#
#
##M1
#sm1<- subset(single.sd, PLASMID=="M1")
##also need to make sure i am only looking at deletions for the amplicon
#m1.del <- subset(s.ap.thirty, PLASMID=="M1")
#toMatch <- substring(sm1$seq, 20) #this removes the first 20 bases in my string to match!
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      m1.del$RECONSTRUCTED_SEQ, value=TRUE)))
##7/15
#matches.all <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                          s.m1$RECONSTRUCTED_SEQ, value=TRUE)))
##14/15? --> but this is for ALL jxns, without a cuttoff. Time to apply the cuttoff then redo it
#s.m1.30 <- subset(s.m1, X._OF_READS>29)
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      s.m1.30$RECONSTRUCTED_SEQ, value=TRUE)))
##12/15
#
##M2
#sm2 <- subset(single.sd, PLASMID=="M2")
##also need to make sure i am only looking at deletions for the amplicon
#m2.del <- subset(s.ap.thirty, PLASMID=="M2")
#toMatch <- sm2$seq
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      m2.del$RECONSTRUCTED_SEQ, value=TRUE)))
##16/19
#matches.all <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                          s.m2$RECONSTRUCTED_SEQ, value=TRUE)))
##20/19? --> but this is for ALL jxns, without a cuttoff. Time to apply the cuttoff then redo it
#s.m2.30 <- subset(s.m2, X._OF_READS>29)
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      s.m2.30$RECONSTRUCTED_SEQ, value=TRUE)))
##17/19! MUCH BETTER! This means that there are several "complex" jxns that are just deletions#
#
##M3
#sm3<- subset(single.sd, PLASMID=="M3")
##also need to make sure i am only looking at deletions for the amplicon
#m3.del <- subset(s.ap.thirty, PLASMID=="M3")
#toMatch <- substring(sm3$seq, 10) #this removes the first 20 bases in my string to match!
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      m3.del$RECONSTRUCTED_SEQ, value=TRUE)))
##9/10
#matches.all <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                          s.m3$RECONSTRUCTED_SEQ, value=TRUE)))
##10/10 --> but this is for ALL jxns, without a cuttoff. Time to apply the cuttoff then redo it
#s.m3.30 <- subset(s.m3, X._OF_READS>29)
#matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
#                                      s.m3.30$RECONSTRUCTED_SEQ, value=TRUE)))
##9/10
#
## >80% of our sanger sequences (Deletions) are found in our amplicon data with a cutoff of 30 reads!

#-----------------------------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------------------------------------------------------
# Changing gears. I need to get some figures to Mitch for the conference.
# 1) distribution of insertions (indels), blunt joins, and mhj

#get all of the insertions
setwd("R:/Varandt/SD-MMEJ/20160616_amplicon_sequence/Combined_Analysis")
a<-read.csv("control_removed_combined_seq-30.csv")
exp.30 <- subset(conr.all.30, SAMPLE=="exp")
exp.30 <- subset(exp.30,CLASS!="exact")
i.30 <- subset(exp.30,PLASMID=="Iw7")
m1.30 <- subset(exp.30,PLASMID=="M1")
m2.30 <- subset(exp.30,PLASMID=="M2")
m3.30 <- subset(exp.30,PLASMID=="M3")

# Need to accurately classify Insertions/Indels. This means that I need to make all of the complex samples that are not just deletions, say insertion, or Indel
toMatch <- iw7.del.all$RECONSTRUCTED_SEQ
i.30$is.del <- regexpr(paste(toMatch,collapse="|"), i.30$RECONSTRUCTED_SEQ)
i.30$CLASS[i.30$is.del==1] <- "deletion"
i.ins.30 <- subset(i.30, is.del==-1)
i.ins.30$is.del <- "InDel"
colnames(i.ins.30)[27] <- "REPAIR_TYPE"
toMatch <- m1.del.all$RECONSTRUCTED_SEQ
m1.30$is.del <- regexpr(paste(toMatch,collapse="|"), m1.30$RECONSTRUCTED_SEQ)
m1.30$CLASS[m1.30$is.del==1] <- "deletion"
m1.ins.30 <- subset(m1.30, is.del==-1)
m1.ins.30$is.del <- "InDel"
colnames(m1.ins.30)[27] <- "REPAIR_TYPE"
toMatch <- m2.del.all$RECONSTRUCTED_SEQ
m2.30$is.del <- regexpr(paste(toMatch,collapse="|"), m2.30$RECONSTRUCTED_SEQ)
m2.30$CLASS[m2.30$is.del==1] <- "deletion"
m2.ins.30 <- subset(m2.30, is.del==-1)
m2.ins.30$is.del <- "InDel"
colnames(m2.ins.30)[27] <- "REPAIR_TYPE"
toMatch <- m3.del.all$RECONSTRUCTED_SEQ
m3.30$is.del <- regexpr(paste(toMatch,collapse="|"), m3.30$RECONSTRUCTED_SEQ)
m3.30$CLASS[m3.30$is.del==1] <- "deletion"
m3.ins.30 <- subset(m3.30, is.del==-1)
m3.ins.30$is.del <- "InDel"
colnames(m3.ins.30)[27] <- "REPAIR_TYPE"

#get SD-MMEJ merged data set to have the columns as the indels
fu <- ip.merge[,2:18]
fu <- cbind(fu, RECONSTRUCTED_SEQ=ip.merge[,1])
fu <- cbind(fu, ip.merge[,19:26])
fu <- cbind(fu, REPAIR_TYPE=ip.merge$REPAIR_TYPE)
names(fu) <- names(i.ins.30)
i.30.RT <- rbind(i.ins.30,fu)
sum(as.numeric(i.30.RT$percent_inaccurate)) # [1] 67.12406% inaccurate reads. This is less than 100, because the rest have been removed due to filtering
fu <- m1p.merge[,2:18]
fu <- cbind(fu, RECONSTRUCTED_SEQ=m1p.merge[,1])
fu <- cbind(fu, m1p.merge[,19:26])
fu <- cbind(fu, REPAIR_TYPE=m1p.merge$REPAIR_TYPE)
names(fu) <- names(m1.ins.30)
m1.30.RT <- rbind(m1.ins.30,fu)
sum(as.numeric(m1.30.RT$percent_inaccurate)) # [1] 59.3564% inaccurate reads. This is less than 100, because the rest have been removed due to filtering
fu <- m2p.merge[,2:18]
fu <- cbind(fu, RECONSTRUCTED_SEQ=m2p.merge[,1])
fu <- cbind(fu, m2p.merge[,19:26])
fu <- cbind(fu, REPAIR_TYPE=m2p.merge$REPAIR_TYPE)
names(fu) <- names(m2.ins.30)
m2.30.RT <- rbind(m2.ins.30,fu)
sum(as.numeric(m2.30.RT$percent_inaccurate)) # [1] 68.50548% inaccurate reads. This is less than 100, because the rest have been removed due to filtering
fu <- m3p.merge[,2:18]
fu <- cbind(fu, RECONSTRUCTED_SEQ=m3p.merge[,1])
fu <- cbind(fu, m3p.merge[,19:26])
fu <- cbind(fu, REPAIR_TYPE=m3p.merge$REPAIR_TYPE)
names(fu) <- names(m3.ins.30)
m3.30.RT <- rbind(m3.ins.30,fu)
sum(as.numeric(m3.30.RT$percent_inaccurate)) # [1] 63.26834% inaccurate reads. This is less than 100, because the rest have been removed due to filtering

# now lets normalize the % of inaccurate reads to the total of reads AFTER filtering
i.30.RT$norm_percent_inaccurate <- i.30.RT$X._OF_READS/sum(i.30.RT$X._OF_READS)*100
sum(i.30.RT$norm_percent_inaccurate)
m1.30.RT$norm_percent_inaccurate <- m1.30.RT$X._OF_READS/sum(m1.30.RT$X._OF_READS)*100
sum(m1.30.RT$norm_percent_inaccurate)
m2.30.RT$norm_percent_inaccurate <- m2.30.RT$X._OF_READS/sum(m2.30.RT$X._OF_READS)*100
sum(m2.30.RT$norm_percent_inaccurate)
m3.30.RT$norm_percent_inaccurate <- m3.30.RT$X._OF_READS/sum(m3.30.RT$X._OF_READS)*100
sum(m3.30.RT$norm_percent_inaccurate)


# merge it all together to make the graph!
all <- rbind(i.30.RT, m1.30.RT, m2.30.RT, m3.30.RT)
sum(all$norm_percent_inaccurate)
write.csv(all, "control_removed_data_percent_inaccurate_repair_type.csv")

#aggregating the data allows me to make better figures
tbls <- aggregate(as.numeric(percent_inaccurate)~REPAIR_TYPE+PLASMID,data=all, sum)
#correct order
x <- subset(tbls, REPAIR_TYPE=="ABJ")
y <- subset(tbls, REPAIR_TYPE=="MHJ")
z <- subset(tbls, REPAIR_TYPE=="InDel")
xyz <- rbind(z,y,x)
colnames(xyz)[3] <- "percent_inaccurate"
levels(xyz$REPAIR_TYPE)
xyz$REPAIR_TYPE <- factor(xyz$REPAIR_TYPE,levels = xyz$REPAIR_TYPE,ordered = TRUE)
type <- ggplot(xyz, aes(fill=REPAIR_TYPE))
type + geom_bar(aes(x=PLASMID, y=percent_inaccurate), position="stack", stat="identity", colour="black") +
  theme_classic(base_size = 18)+
  guides(fill=FALSE)+
  scale_fill_grey(start=0.0, end=1)+
  scale_y_continuous(name = "% Inaccurate Reads",
                     labels=c("0","25","50","75"),
                     expand=c(0,0),
                     limits=c(0,76),
                     breaks=c(0,25,50,75))+
  scale_x_discrete(name = "Plasmid")+
  theme(axis.text.x=element_text(size=18, face="bold"),
        axis.text.y=element_text(size=18, face="bold"),
        axis.title.x=element_text(size=20, face="bold"),
        axis.title.y=element_text(size=20, face="bold"))
#SAVE IT!
# Remake graph but with NORMALIZED Percent Inaccurate Reads
#aggregating the data allows me to make better figures
tbls2 <- aggregate(norm_percent_inaccurate~REPAIR_TYPE+PLASMID,data=all, sum)
#correct order
x2 <- subset(tbls2, REPAIR_TYPE=="ABJ")
y2 <- subset(tbls2, REPAIR_TYPE=="MHJ")
z2 <- subset(tbls2, REPAIR_TYPE=="InDel")
xyz2 <- rbind(z2,y2,x2)
levels(xyz2$REPAIR_TYPE)
xyz2$REPAIR_TYPE <- factor(xyz2$REPAIR_TYPE,levels = xyz2$REPAIR_TYPE,ordered = TRUE)
type <- ggplot(xyz2, aes(fill=REPAIR_TYPE))
type + geom_bar(aes(x=PLASMID, y=norm_percent_inaccurate), position="stack", stat="identity", colour="black") +
  theme_classic(base_size = 18)+
  guides(fill=FALSE)+
  scale_fill_grey(start=0.0, end=1)+
  scale_y_continuous(name = "% Inaccurate Reads",
                     labels=c("0","25","50","75","100"),
                     expand=c(0,0),
                     limits=c(0,105),
                     breaks=c(0,25,50,75,100))+
  scale_x_discrete(name = "Plasmid")+
  theme(axis.text.x=element_text(size=18, face="bold"),
        axis.text.y=element_text(size=18, face="bold"),
        axis.title.x=element_text(size=20, face="bold"),
        axis.title.y=element_text(size=20, face="bold"))

all <- read.csv("control_removed_data_percent_inaccurate_repair_type.csv")
del <- subset(all, re)


#------Insertion Analysis
#Need to get the insertion sequences to run through my R program "INSERTION_PRORGAM_20160725"
# in theory you could just tell the program to pick out this collumn....but i swear to god if you try and do that and it fucks up, dont ask me for help.
# get RECONSTRUCTED_SEQ to run through python SD-MMEJ program
i.ins.seq <- as.data.frame(i.ins.30$RECONSTRUCTED_SEQ)
names(i.ins.seq) <- "RECONSTRUCTED_SEQ"
m1.ins.seq <- as.data.frame(m1.ins.30$RECONSTRUCTED_SEQ)
names(m1.ins.seq) <- "RECONSTRUCTED_SEQ"
m2.ins.seq <- as.data.frame(m2.ins.30$RECONSTRUCTED_SEQ)
names(m2.ins.seq) <- "RECONSTRUCTED_SEQ"
m3.ins.seq <- as.data.frame(m3.ins.30$RECONSTRUCTED_SEQ)
names(m3.ins.seq) <- "RECONSTRUCTED_SEQ"

# save to run through R insertion analysis program. This is a seperate script designed to only need to change certain parameters at the top and then run through. Make sure you have an version that you never go to unless you fuck up something. I dont know how i managed to make this script in the first place, and i dont know how i will be able to fix it.
write.csv(i.ins.seq, "iw7_ins_for_R_script.csv",row.names=FALSE)
write.csv(m1.ins.seq, "m1_ins_for_R_script.csv",row.names=FALSE)
write.csv(m2.ins.seq, "m2_ins_for_R_script.csv",row.names=FALSE)
write.csv(m3.ins.seq, "m3_ins_for_R_script.csv",row.names=FALSE)

#-----------------------------------------------------------------------------
# Ran the through the insertion program for each plasmid. I had to go through manually to change some indexes for the start and end of repeat motifs. Also to change some consistency. Then i reran through the end of the insertion program to remake the tables. They are labeled with a "2" at the end. This step are in the individual R scripts saved in the insertion_analysis folder

# I made an R script that does analysis of insertions that are single step insertions, and determine if they are SD-MMEJ consistent
# Now add these data to the deletion data and do some fucking stats!
setwd("R:/Varandt/SD-MMEJ/20160616_amplicon_sequence/Fixed_MHanalysis_11_20_16")
m <- read.csv("control_removed_data_percent_inaccurate_repair_type.csv")
setwd("R:/Varandt/SD-MMEJ/20160616_amplicon_sequence/Fixed_MHanalysis_11_20_16/insertion_analysis")
a <- read.csv("iw7_insertion_consistency2.csv")
b <- read.csv("m1_insertion_consistency2.csv")
c <- read.csv("m2_insertion_consistency2.csv")
d <- read.csv("m3_insertion_consistency2.csv")

# The deletion boundaries are totally off on the insertions, I need to take into account all of the insertions and shit.
a$RIGHT_DEL <- 165-(nchar(as.character(a$RECONSTRUCTED_SEQ))-a$right_del+1)
b$RIGHT_DEL <- 165-(nchar(as.character(b$RECONSTRUCTED_SEQ))-b$right_del+1)
c$RIGHT_DEL <- 165-(nchar(as.character(c$RECONSTRUCTED_SEQ))-c$right_del+1)
d$RIGHT_DEL <- 165-(nchar(as.character(d$RECONSTRUCTED_SEQ))-d$right_del+1)


iw7 <- subset(m, PLASMID=="Iw7")
m1 <- subset(m, PLASMID=="M1")
m2 <- subset(m, PLASMID=="M2")
m3 <- subset(m, PLASMID=="M3")
iw7ins <- subset(iw7, REPAIR_TYPE=="InDel")
m1ins <- subset(m1, REPAIR_TYPE=="InDel")
m2ins <- subset(m2, REPAIR_TYPE=="InDel")
m3ins <- subset(m3, REPAIR_TYPE=="InDel")

a$plasmid <- "Iw7"
b$plasmid <- "M1"
c$plasmid <- "M2"
d$plasmid <- "M3"
# Need to remove the insertion sequences where it says TRUE and FALSE...they just need to all say TRUE, remove the row that says FALSE
asort <-a[order(a$consistency, decreasing=TRUE),]
ashort <- asort[!duplicated(asort[,8]),]
ashort2 <- ashort[,7:9]
ashort2$right_del <- ashort[,18]
iw7a <- merge(iw7ins, ashort2, by="RECONSTRUCTED_SEQ")
# Calculate percent consistency
max.iip <- sum(iw7a$percent_inaccurate)
max.iip
# total percentage of innacurate reads are insertions or indels = 24.60654%
ins.c <- sum(iw7a[which(iw7a[,30]=="TRUE"),26])
ins.c # 12.51706% of insertion reads are consistent
ins.c/max.iip*100
#[1] 50.86882
# Percent of one step SD-MMEJ consistent insertions is 50.86882%

# M1
# Need to remove the insertion sequences where it says TRUE and FALSE...they just need to all say TRUE, remove the row that says FALSE
bsort <-b[order(b$consistency, decreasing=TRUE),]
bshort <- bsort[!duplicated(bsort[,8]),]
bshort2 <- bshort[,7:9]
bshort2$right_del <- bshort[,18]
m1a <- merge(m1ins, bshort2, by="RECONSTRUCTED_SEQ")
# Calculate percent consistency
max.m1ip <- sum(m1a$percent_inaccurate)
max.m1ip
# total percentage of innacurate reads are insertions or indels = 29.27098%
ins.c <- sum(m1a[which(m1a[,30]=="TRUE"),26])
ins.c # 10.9637
ins.c/max.m1ip*100
#[1] 41.12334
# Percent of one step SD-MMEJ consistent insertions is 37.45586%

# M2
# Need to remove the insertion sequences where it says TRUE and FALSE...they just need to all say TRUE, remove the row that says FALSE
csort <-c[order(c$consistency, decreasing=TRUE),]
cshort <- csort[!duplicated(csort[,8]),]
cshort2 <- cshort[,7:9]
cshort2$right_del <- cshort[,18]
m2a <- merge(m2ins, cshort2, by="RECONSTRUCTED_SEQ")
# Calculate percent consistency
max.m2ip <- sum(m2a$percent_inaccurate)
max.m2ip
# total percentage of innacurate reads are insertions or indels = 19.69838%
ins.c <- sum(m2a[which(m2a[,30]=="TRUE"),26])
ins.c # 10.69846
ins.c/max.m2ip*100
#[1] 54.31136
# Percent of one step SD-MMEJ consistent insertions is 54.89345%

# M3
# Need to remove the insertion sequences where it says TRUE and FALSE...they just need to all say TRUE, remove the row that says FALSE
dsort <-d[order(d$consistency, decreasing=TRUE),]
dshort <- dsort[!duplicated(dsort[,8]),]
dshort2 <- dshort[,7:9]
dshort2$right_del <- dshort[,18]
m3a <- merge(m3ins, dshort2, by="RECONSTRUCTED_SEQ")
# Calculate percent consistency
max.m3ip <- sum(m3a$percent_inaccurate)
max.m3ip
# total percentage of innacurate reads are insertions or indels = 20.95105%
ins.c <- sum(m3a[which(m3a[,30]=="TRUE"),26])
ins.c # 11.52673
ins.c/max.m3ip*100
#[1] 55.01741
# Percent of one step SD-MMEJ consistent insertions is 55.01741%

# Do it for normalized
#Iw7
# Calculate percent consistency
max.iip <- sum(iw7a$norm_percent_inaccurate)
max.iip
# total percentage of innacurate reads are insertions or indels = 36.65831%
ins.c <- sum(iw7a[which(iw7a[,30]=="TRUE"),29])
ins.c # 18.64765% of insertion reads are consistent
ins.c/max.iip*100
#[1] 50.86882
# Percent of one step SD-MMEJ consistent insertions is 50.86882%

#M1
# Calculate percent consistency
max.m1ip <- sum(m1a$norm_percent_inaccurate)
max.m1ip
# total percentage of innacurate reads are insertions or indels = 49.31394%
ins.c <- sum(m1a[which(m1a[,30]=="TRUE"),29])
ins.c # 18.47096
ins.c/max.m1ip*100
#[1] 37.45586
# Percent of one step SD-MMEJ consistent insertions is 37.45586%

#M2
# Calculate percent consistency
max.m2ip <- sum(m2a$norm_percent_inaccurate)
max.m2ip
# total percentage of innacurate reads are insertions or indels = 28.75445%
ins.c <- sum(m2a[which(m2a[,30]=="TRUE"),29])
ins.c # 15.61693
ins.c/max.m2ip*100
#[1]54.31136
# Percent of one step SD-MMEJ consistent insertions is 54.31136%

#M3
# Calculate percent consistency
max.m3ip <- sum(m3a$norm_percent_inaccurate)
max.m3ip
# total percentage of innacurate reads are insertions or indels = 33.1146%
ins.c <- sum(m3a[which(m3a[,30]=="TRUE"),29])
ins.c # 18.21879
ins.c/max.m3ip*100
#[1] 55.01741







#---------------------TEMPLATE PLOT------------------------
# Want to make a cool graph of the repair motifs
setwd("R:/Varandt/SD-MMEJ/20160616_amplicon_sequence/insertion_program")
a <- read.csv("iw7_insertion_consistency2.csv")
setwd("R:/Varandt/SD-MMEJ/20160616_amplicon_sequence/Combined_Analysis")
m <- read.csv("control_removed_data_percent_inaccurate.csv")

iw7 <- subset(m, PLASMID=="Iw7")
iw7ins <- subset(iw7, CLASS=="InDel")

a$plasmid <- "Iw7"

# Need to remove the insertion sequences where it says TRUE and FALSE...they just need to all say TRUE, remove the row that says FALSE
asort <-a[order(a$consistency, decreasing=TRUE),]
ashort <- asort[!duplicated(asort[,8]),]
ashort2 <- ashort[,7:9]
ashort2$right_del <- ashort[,18]
iw7a <- merge(iw7ins, ashort2, by="RECONSTRUCTED_SEQ")


a.t <- subset(a, consistency==TRUE)
iw7a <- merge(iw7ins, a.t, by="RECONSTRUCTED_SEQ")
#need to seperate direct repeat, and reverse complement
dr <- iw7a[,30:32]
dr$percent_inaccurate <- iw7a[,26]
dr <- cbind(dr, iw7a[,36:37])
dr$insertion_length <- nchar(as.character(iw7a[,39]))
dr$RECONSTRUCTED_SEQ <- iw7a$RECONSTRUCTED_SEQ
dr.u <- na.omit(unique(dr))
# This gives CORRECT indexes for templates
dr.u$DR_START_TRUE <- ifelse(dr.u$DR_START>dr.u$right_del, 
                             (165-(nchar(as.character(dr.u$RECONSTRUCTED_SEQ))- dr.u$DR_START-dr.u$insertion_length)),
                             dr.u$DR_START)
dr.u$DR_END_TRUE <- ifelse(dr.u$DR_START>dr.u$right_del, 
                           (165-(nchar(as.character(dr.u$RECONSTRUCTED_SEQ))- dr.u$DR_END-dr.u$insertion_length)),
                           dr.u$DR_END)

# see if I can add the jxns that have the same motif indexes
dr.u$motif <- paste(dr.u$DR_START_TRUE, dr.u$DR_END_TRUE, sep="-")
dr.u.ag <- aggregate(percent_inaccurate~motif,data=dr.u, sum)
dr.u.t <- as.data.frame(table(dr.u$motif))
dr.u.agt <- merge(dr.u.ag, dr.u.t, by.x="motif", by.y="Var1")
dr.um <- merge(dr.u, dr.u.agt, by="motif")
#order the data better
dr.uml <- subset(dr.um, DR_START_TRUE<76) 
dr.umr <- subset(dr.um, DR_START_TRUE>=76) 
dr.uml <- dr.uml[order(dr.uml$DR_START_TRUE),]
dr.umr <- dr.umr[order(dr.umr$DR_START_TRUE, decreasing=TRUE),]
dr.uml$order <- paste(10:(9+nrow(dr.uml)), "-", sep="")   
dr.umr$order <- paste(10:(9+nrow(dr.umr)), "-", sep="")   # if the number of rows is greater than the number of letters, then dont put [1:nrow(...)]
dr.um2 <- rbind(dr.uml, dr.umr)
dr.um2 <- dr.um2[order(dr.um2$order,decreasing=TRUE),]
dr.um3 <- dr.um2[!duplicated(dr.um2$motif),]
y.value <- strsplit("AATTCGGTACATTACCCTGttatCCCTAGCGGCCGCATAGGCCACTAGTGGATCTGGATCC", "")
y.value

neworder <- dr.um3$motif
dr.um4 <- arrange(transform(dr.um3,
                            motif=factor(motif,levels=motif)),motif)
p <- ggplot(dr.um4, aes(x=motif))
p + geom_boxplot(aes(ymin = DR_START_TRUE, lower = DR_START_TRUE, middle = DR_START_TRUE, 
                     upper = DR_END_TRUE, ymax = DR_END_TRUE,
                     colour=percent_inaccurate.y, fill=percent_inaccurate.y), 
                 stat = "identity") + 
  coord_flip()+
  scale_colour_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                         values=c(1.0,0.8,0.6,0.4,0.2,0)) +
  scale_fill_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                       values=c(1.0,0.8,0.6,0.4,0.2,0))+
  scale_y_continuous(limits=c(51,117),
                     breaks=51:116,
                     labels = c("A", "G", "C", "C", "G","A", "A", "T", "T", "C", "G", "G", "T", "A", "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "t", "t", "a", "t", "C", "C", "C", "T", "A", "G", "C", "G", "G", "C", "C", "G", "C", "A", "T", "A", "G", "G", "C", "C", "A", "C", "T", "A", "G", "T", "G", "G", "A", "T", "C", "T", "G", "G", "A", "T", "C", "C"),
                     name = "Nucleotide",
                     expand=c(0,0))+
  theme_classic(base_size = 18)+
  theme(axis.text.y= element_blank())+
  geom_hline(yintercept = 75, colour="red")+
  geom_hline(yintercept = 78,colour="red")

# Now do reverse complement
#need to seperate reverse complement
rc <- as.data.frame(iw7a[,30])
rc <- cbind(rc, iw7a[,33:34])
rc$percent_inaccurate <- iw7a[,26]
rc <- cbind(rc, iw7a[,36:37])
rc$insertion_length <- nchar(as.character(iw7a[,39]))
rc$RECONSTRUCTED_SEQ <- iw7a$RECONSTRUCTED_SEQ
rc.u <- na.omit(unique(rc))
# This gives CORRECT indexes for templates
rc.u$RC_START_TRUE <- ifelse(rc.u$RC_END>rc.u$right_del, 
                             (165-(nchar(as.character(rc.u$RECONSTRUCTED_SEQ))- rc.u$RC_START-rc.u$insertion_length)),
                             rc.u$RC_START)
rc.u$RC_END_TRUE <- ifelse(rc.u$RC_END>rc.u$left_del, 
                           (165-(nchar(as.character(rc.u$RECONSTRUCTED_SEQ))- rc.u$RC_END-rc.u$insertion_length)),
                           rc.u$RC_END)
rc.u$motif <- paste(rc.u$RC_START_TRUE, rc.u$RC_END_TRUE, sep="-")
rc.u.ag <- aggregate(percent_inaccurate~motif,data=rc.u, sum)
rc.u.t <- as.data.frame(table(rc.u$motif))
rc.u.agt <- merge(rc.u.ag, rc.u.t, by.x="motif", by.y="Var1")
rc.um <- merge(rc.u, rc.u.agt, by="motif")
#order the data better
rc.uml <- subset(rc.um, RC_START_TRUE<76) 
rc.umr <- subset(rc.um, RC_START_TRUE>=76) 
rc.uml <- rc.uml[order(rc.uml$RC_START_TRUE),]
rc.umr <- rc.umr[order(rc.umr$RC_START_TRUE, decreasing=TRUE),]
rc.uml$order <- paste(10:(9+nrow(rc.uml)), "-", sep="")   
rc.umr$order <- paste(10:(9+nrow(rc.umr)), "-", sep="")   # if the number of rows is greater than the number of letters, then dont put [1:nrow(...)]
rc.um2 <- rbind(rc.uml, rc.umr)
rc.um2 <- rc.um2[order(rc.um2$order,decreasing=TRUE),]
rc.um3 <- rc.um2[!duplicated(rc.um2$motif),]
y.value <- strsplit("GCCGAATTCGGTACATTACCCTGttatCCCTAGCGGCCGCATAGGCCACTAGTGGATCTGGATCC", "")
y.value

neworder <- rc.um3$motif
rc.um4 <- arrange(transform(rc.um3,
                            motif=factor(motif,levels=motif)),motif)
p <- ggplot(rc.um4, aes(x=motif))
p + geom_boxplot(aes(ymin = RC_START_TRUE, lower = RC_START_TRUE, middle = RC_START_TRUE, 
                     upper = RC_END_TRUE, ymax = RC_END_TRUE,
                     colour=percent_inaccurate.y, fill=percent_inaccurate.y), 
                 stat = "identity") + 
  coord_flip()+
  scale_colour_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                         values=c(1.0,0.8,0.6,0.4,0.2,0)) +
  scale_fill_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                       values=c(1.0,0.8,0.6,0.4,0.2,0))+
  scale_y_continuous(limits=c(51,117),
                     breaks=51:116,
                     labels = c("A", "G", "C", "C", "G","A", "A", "T", "T", "C", "G", "G", "T", "A", "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "t", "t", "a", "t", "C", "C", "C", "T", "A", "G", "C", "G", "G", "C", "C", "G", "C", "A", "T", "A", "G", "G", "C", "C", "A", "C", "T", "A", "G", "T", "G", "G", "A", "T", "C", "T", "G", "G", "A", "T", "C", "C"),
                     name = "Nucleotide",
                     expand=c(0,0))+
  theme_classic(base_size = 18)+
  theme(axis.text.y= element_blank())+
  geom_hline(yintercept = 75, colour="red")+
  geom_hline(yintercept = 78,colour="red")

# Can I combine the graphs and do a facet grid/wrap?
dr.um4$mechanism <- "Loop-out"
rc.um4$mechanism <- "Snap-back"
names(dr.um4) <- names(rc.um4)
comb <- rbind(dr.um4, rc.um4)

p <- ggplot(comb, aes(x=motif))
p + geom_boxplot(aes(ymin = RC_START_TRUE, lower = RC_START_TRUE, middle = RC_START_TRUE, 
                     upper = RC_END_TRUE, ymax = RC_END_TRUE,
                     colour=percent_inaccurate.y, fill=percent_inaccurate.y), 
                 stat = "identity") + 
  coord_flip()+
  facet_grid(~mechanism, scales="free_x")+
  scale_colour_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                         values=c(1.0,0.8,0.6,0.4,0.2,0),
                         guide=FALSE) +
  scale_fill_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                       values=c(1.0,0.8,0.6,0.4,0.2,0),
                       guide_legend(title="% Inaccurate\nReads"))+
  scale_y_continuous(limits=c(51,117),
                     breaks=51:116,
                     labels = c("A", "G", "C", "C", "G","A", "A", "T", "T", "C", "G", "G", "T", "A", "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "t", "t", "a", "t", "C", "C", "C", "T", "A", "G", "C", "G", "G", "C", "C", "G", "C", "A", "T", "A", "G", "G", "C", "C", "A", "C", "T", "A", "G", "T", "G", "G", "A", "T", "C", "T", "G", "G", "A", "T", "C", "C"),
                     name = "Nucleotide",
                     expand=c(0,0))+
  theme_classic(base_size = 18)+
  theme(axis.text.y= element_blank())+
  geom_hline(yintercept = 75, colour="red")+
  geom_hline(yintercept = 78,colour="red")
# not perfect, but save it, to show mitch and terrence

# Try without facet grid
combl <- subset(comb, RC_START_TRUE<76) 
combr <- subset(comb, RC_START_TRUE>=76) 
combl <- combl[order(combl$RC_START_TRUE),]
combr <- combr[order(combr$RC_START_TRUE, decreasing=TRUE),]
combl$order <- paste(10:(9+nrow(combl)), "-", sep="")   
combr$order <- paste(10:(9+nrow(combr)), "-", sep="")   

comb2 <- rbind(combl, combr)
comb2 <- comb2[order(comb2$order,decreasing=TRUE),]
comb3 <- comb2[!duplicated(comb2$motif),]

neworder <- comb3$motif
comb4 <- arrange(transform(comb3,
                           motif=factor(motif,levels=motif)),motif)

p <- ggplot(comb4, aes(x=motif))
p + geom_boxplot(aes(ymin = RC_START_TRUE, lower = RC_START_TRUE, middle = RC_START_TRUE, 
                     upper = RC_END_TRUE, ymax = RC_END_TRUE,
                     fill=percent_inaccurate.y, linetype=mechanism),
                 colour="white",
                 stat = "identity",
                 width=.8) + 
  coord_flip()+
  scale_fill_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                       values=c(1.0,0.8,0.6,0.4,0.2,0),
                       guide_legend(title="% Inaccurate\nReads"))+
  scale_y_continuous(limits=c(51,117),
                     breaks=51:116,
                     labels = c("A", "G", "C", "C", "G","A", "A", "T", "T", "C", "G", "G", "T", "A", "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "t", "t", "a", "t", "C", "C", "C", "T", "A", "G", "C", "G", "G", "C", "C", "G", "C", "A", "T", "A", "G", "G", "C", "C", "A", "C", "T", "A", "G", "T", "G", "G", "A", "T", "C", "T", "G", "G", "A", "T", "C", "C"),
                     name = "Nucleotide",
                     expand=c(0,0))+
  theme_classic(base_size = 18)+
  theme(axis.text.y= element_blank())+
  geom_hline(yintercept = 75, colour="red")+
  geom_hline(yintercept = 78,colour="red")
# save it! "iw7_ins_template-v2.pdf"
write.csv(comb4, "iw7_insertions_template_index.csv")

# TEMPLATE PLOT----------------------------------------------------------------------------------------------------------------------------------
# Mostly good, but i want to switch gears and focus on some insertion stuff. Try and make a graph with all of the insertions for the different plasmids
setwd("R:/Varandt/SD-MMEJ/20160616_amplicon_sequence/insertion_program")
b <- read.csv("m1_insertion_consistency2.csv")
setwd("R:/Varandt/SD-MMEJ/20160616_amplicon_sequence/Combined_Analysis")
m <- read.csv("control_removed_data_percent_inaccurate.csv")

m1 <- subset(m, PLASMID=="M1")
m1ins <- subset(m1, CLASS=="InDel")

a$plasmid <- "m1"

# Need to remove the insertion sequences where it says TRUE and FALSE...they just need to all say TRUE, remove the row that says FALSE
a.t <- subset(a, consistency==TRUE)
m1a <- merge(m1ins, a.t, by="RECONSTRUCTED_SEQ")
#need to seperate direct repeat, and reverse complement
dr <- m1a[,30:32]
dr$percent_inaccurate <- m1a[,26]
dr <- cbind(dr, m1a[,36:37])
dr$insertion_length <- nchar(as.character(m1a[,39]))
dr$RECONSTRUCTED_SEQ <- m1a$RECONSTRUCTED_SEQ
dr.u <- na.omit(unique(dr))
# This gives CORRECT indexes for templates
dr.u$DR_START_TRUE <- ifelse(dr.u$DR_START>dr.u$right_del, 
                             (165-(nchar(as.character(dr.u$RECONSTRUCTED_SEQ))- dr.u$DR_START-dr.u$insertion_length)),
                             dr.u$DR_START)
dr.u$DR_END_TRUE <- ifelse(dr.u$DR_START>dr.u$right_del, 
                           (165-(nchar(as.character(dr.u$RECONSTRUCTED_SEQ))- dr.u$DR_END-dr.u$insertion_length)),
                           dr.u$DR_END)
# see if I can add the jxns that have the same motif indexes
dr.u$motif <- paste(dr.u$DR_START_TRUE, dr.u$DR_END_TRUE, sep="-")
dr.u.ag <- aggregate(percent_inaccurate~motif,data=dr.u, sum)
dr.u.t <- as.data.frame(table(dr.u$motif))
dr.u.agt <- merge(dr.u.ag, dr.u.t, by.x="motif", by.y="Var1")
dr.um <- merge(dr.u, dr.u.agt, by="motif")
#order the data better
dr.uml <- subset(dr.um, DR_START_TRUE<76) 
dr.umr <- subset(dr.um, DR_START_TRUE>=76) 
dr.uml <- dr.uml[order(dr.uml$DR_START_TRUE),]
dr.umr <- dr.umr[order(dr.umr$DR_START_TRUE, decreasing=TRUE),]
dr.uml$order <- paste(10:(9+nrow(dr.uml)), "-", sep="")   
dr.umr$order <- paste(10:(9+nrow(dr.umr)), "-", sep="")   # if the number of rows is greater than the number of letters, then dont put [1:nrow(...)]
dr.um2 <- rbind(dr.uml, dr.umr)
dr.um2 <- dr.um2[order(dr.um2$order,decreasing=TRUE),]
dr.um3 <- dr.um2[!duplicated(dr.um2$motif),]
y.value <- strsplit("AATTCGGTACATTACCCTGttatCCCTAGCGGCCGCATAGGCCACTAGTGGATCTGGATCC", "")
y.value

neworder <- dr.um3$motif
dr.um4 <- arrange(transform(dr.um3,
                            motif=factor(motif,levels=motif)),motif)
p <- ggplot(dr.um4, aes(x=motif))
p + geom_boxplot(aes(ymin = DR_START_TRUE, lower = DR_START_TRUE, middle = DR_START_TRUE, 
                     upper = DR_END_TRUE, ymax = DR_END_TRUE,
                     colour=percent_inaccurate.y, fill=percent_inaccurate.y), 
                 stat = "identity") + 
  coord_flip()+
  scale_colour_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                         values=c(1.0,0.8,0.6,0.4,0.2,0)) +
  scale_fill_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                       values=c(1.0,0.8,0.6,0.4,0.2,0))+
  scale_y_continuous(limits=c(51,117),
                     breaks=51:116,
                     labels = c("A", "G", "C", "C", "G","A", "A", "T", "T", "C", "G", "G", "T", "A", "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "t", "t", "a", "t", "C", "C", "C", "T", "A", "G", "C", "G", "G", "C", "C", "G", "C", "A", "T", "A", "G", "G", "C", "C", "A", "C", "T", "A", "G", "T", "G", "G", "A", "T", "C", "T", "G", "G", "A", "T", "C", "C"),
                     name = "Nucleotide",
                     expand=c(0,0))+
  theme_classic(base_size = 18)+
  theme(axis.text.y= element_blank())+
  geom_hline(yintercept = 75, colour="red")+
  geom_hline(yintercept = 78,colour="red")
# Reverse complement
rc <- as.data.frame(m1a[,30])
rc <- cbind(rc, m1a[,33:34])
rc$percent_inaccurate <- m1a[,26]
rc <- cbind(rc, m1a[,36:37])
rc$insertion_length <- nchar(as.character(m1a[,39]))
rc$RECONSTRUCTED_SEQ <- m1a$RECONSTRUCTED_SEQ
rc.u <- na.omit(unique(rc))
# This gives CORRECT indexes for templates
rc.u$RC_START_TRUE <- ifelse(rc.u$RC_END>rc.u$right_del, 
                             (165-(nchar(as.character(rc.u$RECONSTRUCTED_SEQ))- rc.u$RC_START-rc.u$insertion_length)),
                             rc.u$RC_START)
rc.u$RC_END_TRUE <- ifelse(rc.u$RC_END>rc.u$left_del, 
                           (165-(nchar(as.character(rc.u$RECONSTRUCTED_SEQ))- rc.u$RC_END-rc.u$insertion_length)),
                           rc.u$RC_END)
rc.u$motif <- paste(rc.u$RC_START_TRUE, rc.u$RC_END_TRUE, sep="-")
rc.u.ag <- aggregate(percent_inaccurate~motif,data=rc.u, sum)
rc.u.t <- as.data.frame(table(rc.u$motif))
rc.u.agt <- merge(rc.u.ag, rc.u.t, by.x="motif", by.y="Var1")
rc.um <- merge(rc.u, rc.u.agt, by="motif")
#order the data better
rc.uml <- subset(rc.um, RC_START_TRUE<76) 
rc.umr <- subset(rc.um, RC_START_TRUE>=76) 
rc.uml <- rc.uml[order(rc.uml$RC_START_TRUE),]
rc.umr <- rc.umr[order(rc.umr$RC_START_TRUE, decreasing=TRUE),]
rc.uml$order <- paste(10:(9+nrow(rc.uml)), "-", sep="")   
rc.umr$order <- paste(10:(9+nrow(rc.umr)), "-", sep="")   # if the number of rows is greater than the number of letters, then dont put [1:nrow(...)]
rc.um2 <- rbind(rc.uml, rc.umr)
rc.um2 <- rc.um2[order(rc.um2$order,decreasing=TRUE),]
rc.um3 <- rc.um2[!duplicated(rc.um2$motif),]
y.value <- strsplit("GCCGAATTCGGTACATTACCCTGttatCCCTAGCGGCCGCATAGGCCACTAGTGGATCTGGATCC", "")
y.value

neworder <- rc.um3$motif
rc.um4 <- arrange(transform(rc.um3,
                            motif=factor(motif,levels=motif)),motif)
p <- ggplot(rc.um4, aes(x=motif))
p + geom_boxplot(aes(ymin = RC_START_TRUE, lower = RC_START_TRUE, middle = RC_START_TRUE, 
                     upper = RC_END_TRUE, ymax = RC_END_TRUE,
                     colour=percent_inaccurate.y, fill=percent_inaccurate.y), 
                 stat = "identity") + 
  coord_flip()+
  scale_colour_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                         values=c(1.0,0.8,0.6,0.4,0.2,0)) +
  scale_fill_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                       values=c(1.0,0.8,0.6,0.4,0.2,0))+
  scale_y_continuous(limits=c(51,117),
                     breaks=51:116,
                     labels = c("A", "G", "C", "C", "G","A", "A", "T", "T", "C", "G", "G", "T", "A", "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "t", "t", "a", "t", "C", "C", "C", "T", "A", "G", "C", "G", "G", "C", "C", "G", "C", "A", "T", "A", "G", "G", "C", "C", "A", "C", "T", "A", "G", "T", "G", "G", "A", "T", "C", "T", "G", "G", "A", "T", "C", "C"),
                     name = "Nucleotide",
                     expand=c(0,0))+
  theme_classic(base_size = 18)+
  theme(axis.text.y= element_blank())+
  geom_hline(yintercept = 75, colour="red")+
  geom_hline(yintercept = 78,colour="red")

# Can I combine the graphs and do a facet grid/wrap?
dr.um4$mechanism <- "Loop-out"
rc.um4$mechanism <- "Snap-back"
names(dr.um4) <- names(rc.um4)
m1comb <- rbind(dr.um4, rc.um4)

p <- ggplot(m1comb, aes(x=motif))
p + geom_boxplot(aes(ymin = RC_START_TRUE, lower = RC_START_TRUE, middle = RC_START_TRUE, 
                     upper = RC_END_TRUE, ymax = RC_END_TRUE,
                     colour=percent_inaccurate.y, fill=percent_inaccurate.y), 
                 stat = "identity") + 
  coord_flip()+
  facet_grid(~mechanism, scales="free_x")+
  scale_colour_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                         values=c(1.0,0.8,0.6,0.4,0.2,0),
                         guide=FALSE) +
  scale_fill_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                       values=c(1.0,0.8,0.6,0.4,0.2,0),
                       guide_legend(title="% Inaccurate\nReads"))+
  scale_y_continuous(limits=c(51,117),
                     breaks=51:116,
                     labels = c("A", "G", "C", "C", "G","A", "A", "T", "T", "C", "G", "G", "T", "A", "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "t", "t", "a", "t", "C", "C", "C", "T", "A", "G", "C", "G", "G", "C", "C", "G", "C", "A", "T", "A", "G", "G", "C", "C", "A", "C", "T", "A", "G", "T", "G", "G", "A", "T", "C", "T", "G", "G", "A", "T", "C", "C"),
                     name = "Nucleotide",
                     expand=c(0,0))+
  theme_classic(base_size = 18)+
  theme(axis.text.y= element_blank())+
  geom_hline(yintercept = 75, colour="red")+
  geom_hline(yintercept = 78,colour="red")
# not perfect, but save it, to show mitch and terrence

# Try without facet grid
m1combl <- subset(m1comb, RC_START_TRUE<76) 
m1combr <- subset(m1comb, RC_START_TRUE>=76) 
m1combl <- m1combl[order(m1combl$RC_START_TRUE),]
m1combr <- m1combr[order(m1combr$RC_START_TRUE, decreasing=TRUE),]
m1combl$order <- paste(10:(9+nrow(m1combl)), "-", sep="")   
m1combr$order <- paste(10:(9+nrow(m1combr)), "-", sep="")   

m1comb2 <- rbind(m1combl, m1combr)
m1comb2 <- m1comb2[order(m1comb2$order,decreasing=TRUE),]
m1comb3 <- m1comb2[!duplicated(m1comb2$motif),]

neworder <- m1comb3$motif
m1comb4 <- arrange(transform(m1comb3,
                             motif=factor(motif,levels=motif)),motif)

p <- ggplot(m1comb4, aes(x=motif))
p + geom_boxplot(aes(ymin = RC_START_TRUE, lower = RC_START_TRUE, middle = RC_START_TRUE, 
                     upper = RC_END_TRUE, ymax = RC_END_TRUE,
                     fill=percent_inaccurate.y, linetype=mechanism),
                 colour="white",
                 stat = "identity",
                 width=.8) + 
  coord_flip()+
  scale_fill_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                       values=c(1.0,0.8,0.6,0.4,0.2,0),
                       guide_legend(title="% Inaccurate\nReads"))+
  scale_y_continuous(limits=c(51,117),
                     breaks=51:116,
                     labels = c("A", "G", "C", "C", "G","A", "A", "T", "T", "C", "G", "G", "T", "A", "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "t", "t", "a", "t", "C", "C", "C", "T", "A", "G", "C", "G", "G", "C", "C", "G", "C", "A", "T", "A", "G", "G", "C", "C", "A", "C", "T", "A", "G", "T", "G", "G", "A", "T", "C", "T", "G", "G", "A", "T", "C", "C"),
                     name = "Nucleotide",
                     expand=c(0,0))+
  theme_classic(base_size = 18)+
  theme(axis.text.y= element_blank())+
  geom_hline(yintercept = 75, colour="red")+
  geom_hline(yintercept = 78,colour="red")
# save it! "iw7_ins_template-v2.pdf"
write.csv(m1comb4, "m1_insertions_template_index.csv")

# M2
setwd("R:/Varandt/SD-MMEJ/20160616_amplicon_sequence/insertion_program")
a <- read.csv("m2_insertion_consistency2.csv")
setwd("R:/Varandt/SD-MMEJ/20160616_amplicon_sequence/Combined_Analysis")
m <- read.csv("control_removed_data_percent_inaccurate.csv")

m2 <- subset(m, PLASMID=="M2")
m2ins <- subset(m2, CLASS=="InDel")

a$plasmid <- "m2"

# Need to remove the insertion sequences where it says TRUE and FALSE...they just need to all say TRUE, remove the row that says FALSE
a.t <- subset(a, consistency==TRUE)
m2a <- merge(m2ins, a.t, by="RECONSTRUCTED_SEQ")
#need to seperate direct repeat, and reverse complement
dr <- m2a[,30:32]
dr$percent_inaccurate <- m2a[,26]
dr <- cbind(dr, m2a[,36:37])
dr$insertion_length <- nchar(as.character(m2a[,39]))
dr$RECONSTRUCTED_SEQ <- m2a$RECONSTRUCTED_SEQ
dr.u <- na.omit(unique(dr))
# This gives CORRECT indexes for templates
dr.u$DR_START_TRUE <- ifelse(dr.u$DR_START>dr.u$right_del, 
                             (165-(nchar(as.character(dr.u$RECONSTRUCTED_SEQ))- dr.u$DR_START-dr.u$insertion_length)),
                             dr.u$DR_START)
dr.u$DR_END_TRUE <- ifelse(dr.u$DR_START>dr.u$right_del, 
                           (165-(nchar(as.character(dr.u$RECONSTRUCTED_SEQ))- dr.u$DR_END-dr.u$insertion_length)),
                           dr.u$DR_END)
# see if I can add the jxns that have the same motif indexes
dr.u$motif <- paste(dr.u$DR_START_TRUE, dr.u$DR_END_TRUE, sep="-")
dr.u.ag <- aggregate(percent_inaccurate~motif,data=dr.u, sum)
dr.u.t <- as.data.frame(table(dr.u$motif))
dr.u.agt <- merge(dr.u.ag, dr.u.t, by.x="motif", by.y="Var1")
dr.um <- merge(dr.u, dr.u.agt, by="motif")
#order the data better
dr.uml <- subset(dr.um, DR_START_TRUE<76) 
dr.umr <- subset(dr.um, DR_START_TRUE>=76) 
dr.uml <- dr.uml[order(dr.uml$DR_START_TRUE),]
dr.umr <- dr.umr[order(dr.umr$DR_START_TRUE, decreasing=TRUE),]
dr.uml$order <- paste(10:(9+nrow(dr.uml)), "-", sep="")   
dr.umr$order <- paste(10:(9+nrow(dr.umr)), "-", sep="")   # if the number of rows is greater than the number of letters, then dont put [1:nrow(...)]
dr.um2 <- rbind(dr.uml, dr.umr)
dr.um2 <- dr.um2[order(dr.um2$order,decreasing=TRUE),]
dr.um3 <- dr.um2[!duplicated(dr.um2$motif),]
neworder <- dr.um3$motif
dr.um4 <- arrange(transform(dr.um3,
                            motif=factor(motif,levels=motif)),motif)
#Reverse complement
rc <- as.data.frame(m2a[,30])
rc <- cbind(rc, m2a[,33:34])
rc$percent_inaccurate <- m2a[,26]
rc <- cbind(rc, m2a[,36:37])
rc$insertion_length <- nchar(as.character(m2a[,39]))
rc$RECONSTRUCTED_SEQ <- m2a$RECONSTRUCTED_SEQ
rc.u <- na.omit(unique(rc))
# This gives CORRECT indexes for templates
rc.u$RC_START_TRUE <- ifelse(rc.u$RC_END>rc.u$right_del, 
                             (165-(nchar(as.character(rc.u$RECONSTRUCTED_SEQ))- rc.u$RC_START-rc.u$insertion_length)),
                             rc.u$RC_START)
rc.u$RC_END_TRUE <- ifelse(rc.u$RC_END>rc.u$left_del, 
                           (165-(nchar(as.character(rc.u$RECONSTRUCTED_SEQ))- rc.u$RC_END-rc.u$insertion_length)),
                           rc.u$RC_END)
rc.u$motif <- paste(rc.u$RC_START_TRUE, rc.u$RC_END_TRUE, sep="-")
rc.u.ag <- aggregate(percent_inaccurate~motif,data=rc.u, sum)
rc.u.t <- as.data.frame(table(rc.u$motif))
rc.u.agt <- merge(rc.u.ag, rc.u.t, by.x="motif", by.y="Var1")
rc.um <- merge(rc.u, rc.u.agt, by="motif")
#order the data better
rc.uml <- subset(rc.um, RC_START_TRUE<76) 
rc.umr <- subset(rc.um, RC_START_TRUE>=76) 
rc.uml <- rc.uml[order(rc.uml$RC_START_TRUE),]
rc.umr <- rc.umr[order(rc.umr$RC_START_TRUE, decreasing=TRUE),]
rc.uml$order <- paste(10:(9+nrow(rc.uml)), "-", sep="")   
rc.umr$order <- paste(10:(9+nrow(rc.umr)), "-", sep="")   # if the number of rows is greater than the number of letters, then dont put [1:nrow(...)]
rc.um2 <- rbind(rc.uml, rc.umr)
rc.um2 <- rc.um2[order(rc.um2$order,decreasing=TRUE),]
rc.um3 <- rc.um2[!duplicated(rc.um2$motif),]
neworder <- rc.um3$motif
rc.um4 <- arrange(transform(rc.um3,
                            motif=factor(motif,levels=motif)),motif)
dr.um4$mechanism <- "Loop-out"
rc.um4$mechanism <- "Snap-back"
names(dr.um4) <- names(rc.um4)
m2comb <- rbind(dr.um4, rc.um4)
m2combl <- subset(m2comb, RC_START_TRUE<76) 
m2combr <- subset(m2comb, RC_START_TRUE>=76) 
m2combl <- m2combl[order(m2combl$RC_START_TRUE),]
m2combr <- m2combr[order(m2combr$RC_START_TRUE, decreasing=TRUE),]
m2combl$order <- paste(10:(9+nrow(m2combl)), "-", sep="")   
m2combr$order <- paste(10:(9+nrow(m2combr)), "-", sep="")   

m2comb2 <- rbind(m2combl, m2combr)
m2comb2 <- m2comb2[order(m2comb2$order,decreasing=TRUE),]
m2comb3 <- m2comb2[!duplicated(m2comb2$motif),]

neworder <- m2comb3$motif
m2comb4 <- arrange(transform(m2comb3,
                             motif=factor(motif,levels=motif)),motif)

# M3
setwd("R:/Varandt/SD-MMEJ/20160616_amplicon_sequence/insertion_program")
a <- read.csv("m3_insertion_consistency2.csv")
setwd("R:/Varandt/SD-MMEJ/20160616_amplicon_sequence/Combined_Analysis")
m <- read.csv("control_removed_data_percent_inaccurate.csv")

m3 <- subset(m, PLASMID=="M3")
m3ins <- subset(m3, CLASS=="InDel")

a$plasmid <- "m3"

# Need to remove the insertion sequences where it says TRUE and FALSE...they just need to all say TRUE, remove the row that says FALSE
a.t <- subset(a, consistency==TRUE)
m3a <- merge(m3ins, a.t, by="RECONSTRUCTED_SEQ")
#need to seperate direct repeat, and reverse complement
dr <- m3a[,30:32]
dr$percent_inaccurate <- m3a[,26]
dr <- cbind(dr, m3a[,36:37])
dr$insertion_length <- nchar(as.character(m3a[,39]))
dr$RECONSTRUCTED_SEQ <- m3a$RECONSTRUCTED_SEQ
dr.u <- na.omit(unique(dr))
# This gives CORRECT indexes for templates
dr.u$DR_START_TRUE <- ifelse(dr.u$DR_START>dr.u$right_del, 
                             (165-(nchar(as.character(dr.u$RECONSTRUCTED_SEQ))- dr.u$DR_START-dr.u$insertion_length)),
                             dr.u$DR_START)
dr.u$DR_END_TRUE <- ifelse(dr.u$DR_START>dr.u$right_del, 
                           (165-(nchar(as.character(dr.u$RECONSTRUCTED_SEQ))- dr.u$DR_END-dr.u$insertion_length)),
                           dr.u$DR_END)
# see if I can add the jxns that have the same motif indexes
dr.u$motif <- paste(dr.u$DR_START_TRUE, dr.u$DR_END_TRUE, sep="-")
dr.u.ag <- aggregate(percent_inaccurate~motif,data=dr.u, sum)
dr.u.t <- as.data.frame(table(dr.u$motif))
dr.u.agt <- merge(dr.u.ag, dr.u.t, by.x="motif", by.y="Var1")
dr.um <- merge(dr.u, dr.u.agt, by="motif")
#order the data better
dr.uml <- subset(dr.um, DR_START_TRUE<76) 
dr.umr <- subset(dr.um, DR_START_TRUE>=76) 
dr.uml <- dr.uml[order(dr.uml$DR_START_TRUE),]
dr.umr <- dr.umr[order(dr.umr$DR_START_TRUE, decreasing=TRUE),]
dr.uml$order <- paste(10:(9+nrow(dr.uml)), "-", sep="")   
dr.umr$order <- paste(10:(9+nrow(dr.umr)), "-", sep="")   # if the number of rows is greater than the number of letters, then dont put [1:nrow(...)]
dr.um2 <- rbind(dr.uml, dr.umr)
dr.um2 <- dr.um2[order(dr.um2$order,decreasing=TRUE),]
dr.um3 <- dr.um2[!duplicated(dr.um2$motif),]
neworder <- dr.um3$motif
dr.um4 <- arrange(transform(dr.um3,
                            motif=factor(motif,levels=motif)),motif)
#Reverse complement
rc <- as.data.frame(m3a[,30])
rc <- cbind(rc, m3a[,33:34])
rc$percent_inaccurate <- m3a[,26]
rc <- cbind(rc, m3a[,36:37])
rc$insertion_length <- nchar(as.character(m3a[,39]))
rc$RECONSTRUCTED_SEQ <- m3a$RECONSTRUCTED_SEQ
rc.u <- na.omit(unique(rc))
# This gives CORRECT indexes for templates
rc.u$RC_START_TRUE <- ifelse(rc.u$RC_END>rc.u$right_del, 
                             (165-(nchar(as.character(rc.u$RECONSTRUCTED_SEQ))- rc.u$RC_START-rc.u$insertion_length)),
                             rc.u$RC_START)
rc.u$RC_END_TRUE <- ifelse(rc.u$RC_END>rc.u$left_del, 
                           (165-(nchar(as.character(rc.u$RECONSTRUCTED_SEQ))- rc.u$RC_END-rc.u$insertion_length)),
                           rc.u$RC_END)
rc.u$motif <- paste(rc.u$RC_START_TRUE, rc.u$RC_END_TRUE, sep="-")
rc.u.ag <- aggregate(percent_inaccurate~motif,data=rc.u, sum)
rc.u.t <- as.data.frame(table(rc.u$motif))
rc.u.agt <- merge(rc.u.ag, rc.u.t, by.x="motif", by.y="Var1")
rc.um <- merge(rc.u, rc.u.agt, by="motif")
#order the data better
rc.uml <- subset(rc.um, RC_START_TRUE<76) 
rc.umr <- subset(rc.um, RC_START_TRUE>=76) 
rc.uml <- rc.uml[order(rc.uml$RC_START_TRUE),]
rc.umr <- rc.umr[order(rc.umr$RC_START_TRUE, decreasing=TRUE),]
rc.uml$order <- paste(10:(9+nrow(rc.uml)), "-", sep="")   
rc.umr$order <- paste(10:(9+nrow(rc.umr)), "-", sep="")   # if the number of rows is greater than the number of letters, then dont put [1:nrow(...)]
rc.um2 <- rbind(rc.uml, rc.umr)
rc.um2 <- rc.um2[order(rc.um2$order,decreasing=TRUE),]
rc.um3 <- rc.um2[!duplicated(rc.um2$motif),]
neworder <- rc.um3$motif
rc.um4 <- arrange(transform(rc.um3,
                            motif=factor(motif,levels=motif)),motif)
dr.um4$mechanism <- "Loop-out"
rc.um4$mechanism <- "Snap-back"
names(dr.um4) <- names(rc.um4)
m3comb <- rbind(dr.um4, rc.um4)
m3combl <- subset(m3comb, RC_START_TRUE<76) 
m3combr <- subset(m3comb, RC_START_TRUE>=76) 
m3combl <- m3combl[order(m3combl$RC_START_TRUE),]
m3combr <- m3combr[order(m3combr$RC_START_TRUE, decreasing=TRUE),]
m3combl$order <- paste(10:(9+nrow(m3combl)), "-", sep="")   
m3combr$order <- paste(10:(9+nrow(m3combr)), "-", sep="")   

m3comb2 <- rbind(m3combl, m3combr)
m3comb2 <- m3comb2[order(m3comb2$order,decreasing=TRUE),]
m3comb3 <- m3comb2[!duplicated(m3comb2$motif),]

neworder <- m3comb3$motif
m3comb4 <- arrange(transform(m3comb3,
                             motif=factor(motif,levels=motif)),motif)





# Combine m1 and Iw7, to see what it looks like
comb4$plasmid <- "iw7"
m1comb4$plasmid <- "m1"
m2comb4$plasmid <- "m2"
m3comb4$plasmid <- "m3"
names(comb4) <- names(m1comb4)
names(m2comb4) <- names(m1comb4)
names(m3comb4) <- names(m1comb4)
v <- rbind(comb4, m1comb4, m2comb4,m3comb4)
vl <- subset(v, RC_START_TRUE<76) 
vr <- subset(v, RC_START_TRUE>=76) 
vl <- vl[order(vl$RC_START_TRUE),]
vr <- vr[order(vr$RC_START_TRUE, decreasing=TRUE),]
vl$order <- paste(10:(9+nrow(vl)), "-", sep="")   
vr$order <- paste(10:(9+nrow(vr)), "-", sep="")   

v2 <- rbind(vl, vr)
v2 <- v2[order(v2$order,decreasing=TRUE),]

neworder <- v2$motif
v3 <- arrange(transform(v2,
                        motif=factor(motif,levels=motif)),motif)
summary(v)
#graph it
p <- ggplot(v, aes(x=motif))
p + geom_boxplot(aes(ymin = RC_START_TRUE, lower = RC_START_TRUE, middle = RC_START_TRUE, 
                     upper = RC_END_TRUE, ymax = RC_END_TRUE,
                     fill=percent_inaccurate.y, linetype=mechanism),
                 colour="white",
                 stat = "identity",
                 width=.8) + 
  coord_flip()+
  facet_wrap(~plasmid, ncol=2)+
  scale_fill_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                       values=c(1.0,0.8,0.6,0.4,0.2,0),
                       guide_legend(title="% Inaccurate\nReads"))+
  scale_y_continuous(limits=c(51,117),
                     breaks=51:116,
                     labels = c("A", "G", "C", "C", "G","A", "A", "T", "T", "C", "G", "G", "T", "A", "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "t", "t", "a", "t", "C", "C", "C", "T", "A", "G", "C", "G", "G", "C", "C", "G", "C", "A", "T", "A", "G", "G", "C", "C", "A", "C", "T", "A", "G", "T", "G", "G", "A", "T", "C", "T", "G", "G", "A", "T", "C", "C"),
                     name = "Nucleotide",
                     expand=c(0,0))+
  theme_classic(base_size = 18)+
  theme(axis.text.y= element_blank())+
  geom_hline(yintercept = 75, colour="red")+
  geom_hline(yintercept = 78,colour="red")

# make the graphs seperately, and place into one, like i did with the deletion stoppers
library(grid)
f1 <- ggplot(comb4, aes(x=motif))+
  geom_boxplot(aes(ymin = RC_START_TRUE, lower = RC_START_TRUE, middle = RC_START_TRUE, 
                   upper = RC_END_TRUE, ymax = RC_END_TRUE,
                   fill=percent_inaccurate.y, linetype=mechanism),
               colour="white",
               stat = "identity",
               width=.8) + 
  coord_flip()+
  facet_wrap(~plasmid, ncol=2)+
  scale_fill_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                       values=c(1.0,0.8,0.6,0.4,0.2,0),
                       guide_legend(title="% Inaccurate\nReads"),
                       limits=c(0.03377,3.27715))+
  scale_y_continuous(limits=c(51,117),
                     breaks=51:116,
                     labels = c("A", "G", "C", "C", "G","A", "A", "T", "T", "C", "G", "G", "T", "A", "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "t", "t", "a", "t", "C", "C", "C", "T", "A", "G", "C", "G", "G", "C", "C", "G", "C", "A", "T", "A", "G", "G", "C", "C", "A", "C", "T", "A", "G", "T", "G", "G", "A", "T", "C", "T", "G", "G", "A", "T", "C", "C"),
                     name = "Nucleotide",
                     expand=c(0,0))+
  theme_classic(base_size = 18)+
  theme(axis.text.y= element_blank(),
        axis.title.y=element_blank(),
        axis.ticks=element_blank(),
        panel.grid.minor=element_line(colour="grey", size=0.5))+
  geom_hline(yintercept = c(75,78), colour="red", size=0.75)
f2 <- ggplot(m1comb4, aes(x=motif))+
  geom_boxplot(aes(ymin = RC_START_TRUE, lower = RC_START_TRUE, middle = RC_START_TRUE, 
                   upper = RC_END_TRUE, ymax = RC_END_TRUE,
                   fill=percent_inaccurate.y, linetype=mechanism),
               colour="white",
               stat = "identity",
               width=.8) + 
  coord_flip()+
  facet_wrap(~plasmid, ncol=2)+
  scale_fill_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                       values=c(1.0,0.8,0.6,0.4,0.2,0),
                       guide_legend(title="% Inaccurate\nReads"),
                       limits=c(0.03377,3.27715))+
  scale_y_continuous(limits=c(51,117),
                     breaks=51:116,
                     labels = c("A", "G", "C", "C", "G","A", "A", "T", "T", "C", "G", "G", "T", "A", "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "t", "t", "a", "t", "C", "C", "C", "T", "A", "G", "C", "G", "G", "A", "C", "G", "C", "A", "T", "A", "G", "G", "C", "C", "A", "C", "T", "A", "G", "T", "G", "G", "A", "T", "C", "T", "G", "G", "A", "T", "C", "C"),
                     name = "Nucleotide",
                     expand=c(0,0))+
  theme_classic(base_size = 18)+
  theme(axis.text.y= element_blank(),
        axis.title.y=element_blank(),
        axis.ticks=element_blank(),
        panel.grid.minor=element_line(colour="grey", size=0.5))+
  geom_hline(yintercept = c(75,78), colour="red", size=0.75)
f3 <- ggplot(m2comb4, aes(x=motif))+
  geom_boxplot(aes(ymin = RC_START_TRUE, lower = RC_START_TRUE, middle = RC_START_TRUE, 
                   upper = RC_END_TRUE, ymax = RC_END_TRUE,
                   fill=percent_inaccurate.y, linetype=mechanism),
               colour="white",
               stat = "identity",
               width=.8) + 
  coord_flip()+
  facet_wrap(~plasmid, ncol=2)+
  scale_fill_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                       values=c(1.0,0.8,0.6,0.4,0.2,0),
                       guide_legend(title="% Inaccurate\nReads"),
                       limits=c(0.03377,3.27715))+
  scale_y_continuous(limits=c(51,117),
                     breaks=51:116,
                     labels = c("A", "G", "C", "C", "G","A", "A", "T", "T", "C", "G", "G", "T", "A", "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "t", "t", "a", "t", "C", "C", "C", "T", "A", "G", "T", "G", "G", "A", "C", "G", "C", "A", "T", "A", "G", "G", "C", "C", "A", "C", "T", "A", "G", "T", "G", "G", "A", "T", "C", "T", "G", "G", "A", "T", "C", "C"),
                     name = "Nucleotide",
                     expand=c(0,0))+
  theme_classic(base_size = 18)+
  theme(axis.text.y= element_blank(),
        axis.title.y=element_blank(),
        axis.ticks=element_blank(),
        panel.grid.minor=element_line(colour="grey", size=0.5))+
  geom_hline(yintercept = c(75,78), colour="red", size=0.75)
f4 <- ggplot(m3comb4, aes(x=motif))+
  geom_boxplot(aes(ymin = RC_START_TRUE, lower = RC_START_TRUE, middle = RC_START_TRUE, 
                   upper = RC_END_TRUE, ymax = RC_END_TRUE,
                   fill=percent_inaccurate.y, linetype=mechanism),
               colour="white",
               stat = "identity",
               width=.8) + 
  coord_flip()+
  facet_wrap(~plasmid, ncol=2)+
  scale_fill_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
                       values=c(1.0,0.8,0.6,0.4,0.2,0),
                       guide_legend(title="% Inaccurate\nReads"),
                       limits=c(0.03377,3.27715))+
  scale_y_continuous(limits=c(51,117),
                     breaks=51:116,
                     labels = c("A", "G", "C", "C", "G","A", "A", "T", "T", "C", "G", "G", "T", "A", "C", "A", "T", "T", "A", "C", "C", "C", "T", "G", "t", "t", "a", "t", "C", "C", "C", "T", "A", "G", "T", "G", "G", "A", "C", "G", "C", "A", "T", "A", "G", "G", "C", "C", "A", "C", "T", "A", "G", "T", "G", "G", "A", "T", "C", "T", "G", "G", "A", "T", "C", "C"),
                     name = "Nucleotide",
                     expand=c(0,0))+
  theme_classic(base_size = 18)+
  theme(axis.text.y= element_blank(),
        axis.title.y=element_blank(),
        axis.ticks=element_blank(),
        panel.grid.minor=element_line(colour="grey", size=0.5))+
  geom_hline(yintercept = c(75,78), colour="red", size=0.75)

# saved all of them seperately. they still need some work to make completely pretty, but not bad

#==================================================================================================================
#=====================================================================================================================
# NICE!
# make same graph for sanger data
# Fortunately, i did all the hard work a while ago and make this data in long format...it doesnt have the actual sequences on it, but i will have to get that later
setwd("C:/Users/Varandt/OneDrive/McVey Lab/Projects/SD-MMEJ/Data/New_data_analysis")
injection_data <- read.csv("injection_data.csv")
inj.i <- subset(injection_data, plasmid=="Iw7")
inj.m1<- subset(injection_data, plasmid=="M1")
inj.m2<- subset(injection_data, plasmid=="M2")
inj.m3<- subset(injection_data, plasmid=="M3")
inj.i$count <- 1
inj.i$percent <- inj.i$count/sum(inj.i$count)*100
inj.m1$count <- 1
inj.m1$percent <- inj.m1$count/sum(inj.m1$count)*100
inj.m2$count <- 1
inj.m2$percent <- inj.m2$count/sum(inj.m2$count)*100
inj.m3$count <- 1
inj.m3$percent <- inj.m3$count/sum(inj.m3$count)*100
inj.all <- rbind(inj.i, inj.m1, inj.m2, inj.m3)
inj.indel <- subset(inj.all, Repair_Type=="Insertion") 
inj.abj <- subset(inj.all, Repair_Type=="ABJ") 
inj.mhj <- subset(inj.all, Repair_Type=="MHJ") 

inj.order <- rbind(inj.indel, inj.mhj, inj.abj)

repair_table <- table(inj.order$Repair_Type)

repair_levels <- names(repair_table)[order(repair_table)]
inj.order$Repair_Type2 <- factor(inj.order$Repair_Type, levels = repair_levels)

type.sanger <- ggplot(inj.order, aes(fill=Repair_Type2))
type.sanger + geom_bar(aes(x=plasmid, y=percent), position="stack", stat="identity") +
  theme_bw()
# SAVE IT!







# make new graph of deletion lengths for Deletion class
# make big table with the sd-mmej consistency for all plasmids
all <- rbind(ip.merge, m1p.merge, m2p.merge, m3p.merge)
m.del <- ggplot(all, aes(fill=CONSISTENCY))
m.del <- m.del+geom_bar(aes(x=TOTAL_DELETION.y, y=percent_inaccurate), position="stack", stat="identity")+facet_grid(PLASMID.x~ .)
m.del
update_labels(m.del, list(x = "Deletion Length", y="Percent of Inaccurate Reads"))

#------------------------------------------------------------------------------

# ran through SD-MMEJ consistency program - then tacked them on to the csv file along with the rest of the iw7 samples
setwd("R:/Varandt/SD-MMEJ/20160616_amplicon_sequence/Combined_Analysis/Python_Analysis")
p <- read.csv("06282016_amplicon_python_consistency.csv")
p$SEQ <- toupper(p$SEQ)
# seperate based on plasmid
ip <- subset(p, PLASMID=="Iw7")
m1p <- subset(p, PLASMID=="M1")
m2p <- subset(p, PLASMID=="M2")
test <- subset(ip, !duplicated(ip$SEQ))

# Now calculate percentages
# I need to merge with the del and complex del data set
iw7.mdel <- rbind(s.iw7.del, del[,2:25])

ip.merge <- merge(iw7.mdel, test, by.x="RECONSTRUCTED_SEQ", by.y="SEQ", all=FALSE)
# YAY it works! Now figure out how to do weighted counts or something
max.ip <- sum(ip.merge$percent)
max.ip












s.ap.del <- subset(s.ap, CLASS=="deletion")
summary(s.ap.del)
#only samples with >30 reads
s.ap.thirty <- subset(s.ap.del, X._OF_READS>29)
tail(s.ap.thirty)
summary(s.ap.thirty)
#plot percent per deletion length
attach(s.ap.thirty)
m.del <- ggplot(s.ap.thirty, aes(x=-TOTAL_DELETION, y=percent))
m.del+geom_point(position="stack")+facet_grid(PLASMID~ .)
m.del <- ggplot(s.ap.thirty, aes(x=-TOTAL_DELETION, y=percent))
m.del+geom_line(position="stack")+facet_grid(PLASMID~ .)

#THIS IS WHAT I WANT! BARPLOT/HISOGRAM based on percentages
m.del <- ggplot(s.ap.thirty, aes(x=-TOTAL_DELETION, y=percent, fill=PLASMID))
m.del <- m.del+geom_bar(position="stack", stat="identity")+facet_grid(PLASMID~ .)
update_labels(m.del, list(x = "Deletion Length", y="Percent of Reads"))
# Saved as "amplicon_total_deletion.pdf"
#SAVE IT!

# NOW! Do the same thing i did before with left and right with percentages...BUT, unlike total deletion, I want the multiple alignments, because we want to know both possibilities
# This means I need to figure out how to calculate percentage of a the "duplicated" list.
# I think i know. Just take the sum of reads from s.plasmid, and divide that by the non duplicated --> this will mean than the percentages will NOT add to 100 in the duplicated list, but that is fine

r.con <- sum(s.con$X._OF_READS)
con$percent <-con$X._OF_READS/r.con*100
head(con)
summary(con)
summary(s.con)
# PERFECT! do it for all of them!
r.iw7 <- sum(s.iw7$X._OF_READS)
iw7$percent <-iw7$X._OF_READS/r.iw7*100
head(iw7)
head(s.iw7)
r.m1 <- sum(s.m1$X._OF_READS)
m1$percent <-m1$X._OF_READS/r.m1*100
head(m1)
head(s.m1)
r.m2 <- sum(s.m2$X._OF_READS)
m2$percent <-m2$X._OF_READS/r.m2*100
head(m2)
head(s.m2)
r.m3 <- sum(s.m3$X._OF_READS)
m3$percent <-m3$X._OF_READS/r.m3*100
head(m3)
head(s.m3)
ap <- rbind(con, iw7, m1, m2, m3)
ap.del <- subset(ap, CLASS=="deletion")
summary(ap.del)
#only samples with >30 reads
ap.thirty <- subset(ap.del, X._OF_READS>29)
tail(ap.thirty)
summary(ap.thirty)
#plot percent per deletion length
attach(ap.thirty)
m.del <- ggplot(ap.thirty, aes(fill=PLASMID))
m.del <- m.del+geom_bar(aes(x=-DELETION_FROM_RIGHT, y=percent), position="stack", stat="identity")+geom_bar(aes(x=DELETION_FROM_LEFT, y=percent), position="stack", stat="identity")+facet_grid(PLASMID~ .)
update_labels(m.del, list(x = "Deletion Stopping Position", y="Percent of Reads"))
#SAVE IT! "amplicon_deletion_stopper.pdf"

#I should probably start to compare these to the sanger data!
setwd("R:/Varandt/SD-MMEJ/20160616_amplicon_sequence")
sd<-read.csv("sanger_del.csv")
summary(sd)
m.del <- ggplot(sd, aes(x=TOTAL_DELETION, fill=PLASMID))
m.del+geom_bar(position="stack")+facet_grid(PLASMID~ .)
# try it slightly differently to make the two plots almost identicle
sd$COUNT <- 1
m.del <- ggplot(sd, aes(x=TOTAL_DELETION, y=COUNT, fill=PLASMID))
m.del <- m.del+geom_bar(position="stack", stat="identity")+facet_grid(PLASMID~ .)
update_labels(m.del, list(x = "Deletion Length", y="Number of Occurrences"))
#This is much better. SAVE!

# Need to determine the % of sanger jxns in our amplicon data
# First create no duplicate sanger data
# to do this, i need to remove the '----' from the aligned sequences that were placed into python
single.sd <- subset(sd,!duplicated(sd$SEQ))
head(single.sd)
#create "seq" without "---"
# to do this, i need to remove the '----' from the aligned sequences that were placed into python
single.sd$seq <- gsub("-", "",single.sd$SEQ)
siw7 <- subset(single.sd, PLASMID=="Iw7")
#also need to make sure i am only looking at deletions for the amplicon
iw7.del <- subset(s.ap.thirty, PLASMID=="Iw7")
toMatch <- siw7$seq
matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
                                      iw7.del$RECONSTRUCTED_SEQ, value=TRUE)))
#14/21
matches.all <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
                                          s.iw7$RECONSTRUCTED_SEQ, value=TRUE)))
#25/21? --> but this is for ALL jxns, without a cuttoff. Time to apply the cuttoff then redo it
s.iw7.30 <- subset(s.iw7, X._OF_READS>29)
matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
                                      s.iw7.30$RECONSTRUCTED_SEQ, value=TRUE)))
#20/21! MUCH BETTER! This means that there are several "complex" jxns that are just deletions

#M1
sm1<- subset(single.sd, PLASMID=="M1")
#also need to make sure i am only looking at deletions for the amplicon
m1.del <- subset(s.ap.thirty, PLASMID=="M1")
toMatch <- substring(sm1$seq, 20) #this removes the first 20 bases in my string to match!
matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
                                      m1.del$RECONSTRUCTED_SEQ, value=TRUE)))
#7/15
matches.all <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
                                          s.m1$RECONSTRUCTED_SEQ, value=TRUE)))
#14/15? --> but this is for ALL jxns, without a cuttoff. Time to apply the cuttoff then redo it
s.m1.30 <- subset(s.m1, X._OF_READS>29)
matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
                                      s.m1.30$RECONSTRUCTED_SEQ, value=TRUE)))
#12/15

#M2
sm2 <- subset(single.sd, PLASMID=="M2")
#also need to make sure i am only looking at deletions for the amplicon
m2.del <- subset(s.ap.thirty, PLASMID=="M2")
toMatch <- sm2$seq
matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
                                      m2.del$RECONSTRUCTED_SEQ, value=TRUE)))
#16/19
matches.all <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
                                          s.m2$RECONSTRUCTED_SEQ, value=TRUE)))
#20/19? --> but this is for ALL jxns, without a cuttoff. Time to apply the cuttoff then redo it
s.m2.30 <- subset(s.m2, X._OF_READS>29)
matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
                                      s.m2.30$RECONSTRUCTED_SEQ, value=TRUE)))
#17/19! MUCH BETTER! This means that there are several "complex" jxns that are just deletions

#M3
sm3<- subset(single.sd, PLASMID=="M3")
#also need to make sure i am only looking at deletions for the amplicon
m3.del <- subset(s.ap.thirty, PLASMID=="M3")
toMatch <- substring(sm3$seq, 10) #this removes the first 20 bases in my string to match!
matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
                                      m3.del$RECONSTRUCTED_SEQ, value=TRUE)))
#9/10
matches.all <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
                                          s.m3$RECONSTRUCTED_SEQ, value=TRUE)))
#10/10 --> but this is for ALL jxns, without a cuttoff. Time to apply the cuttoff then redo it
s.m3.30 <- subset(s.m3, X._OF_READS>29)
matches <- as.data.frame(unique (grep(paste(toMatch,collapse="|"), 
                                      s.m3.30$RECONSTRUCTED_SEQ, value=TRUE)))
#9/10

# >80% of our sanger sequences (Deletions) are found in our amplicon data with a cutoff of 30 reads!

#Now I need to get the deletion sequences to put through the python program
# Need to make sure i have version with NO DUPLICATES!

i.seq <- as.data.frame(iw7.del$ALIGNED_SEQ)
m1.seq <- as.data.frame(m1.del$ALIGNED_SEQ)
m2.seq <- as.data.frame(m2.del$ALIGNED_SEQ)
m3.seq <- as.data.frame(m3.del$ALIGNED_SEQ) # This one doesnt actually have the aligned output for osme reason

#Now save them
setwd("R:/Varandt/SD-MMEJ/20160616_amplicon_sequence/Combined_Analysis")
write.csv(i.seq, file="iw7_del_seq_for_python.csv")
write.csv(m1.seq, file="m1_del_seq_for_python.csv")
write.csv(m2.seq, file="m2_del_seq_for_python.csv")
write.csv(m3.seq, file="m3_del_seq_for_python.csv")

#Ran through sequences through the SD-MMEJ python program and got back the data.
# I have raw counts for consistency, but i need to correct them based on percent of reads.
# from the output file, i made a new table with that has jxn and consistency to merge with the input file

setwd("R:/Varandt/SD-MMEJ/20160616_amplicon_sequence/Combined_Analysis/Python_Analysis")
ip <- read.csv("06272016_iw7_amplicon_merge.csv")
#need to make all characters upper case
ip$SEQ <- toupper(ip$SEQ) 
ip.merge <- merge(iw7.del, ip, by.x="RECONSTRUCTED_SEQ", by.y="SEQ", all=FALSE)
# for some reason there are two that were duplicated...dont know why
s.exp2 <- subset(exp, !duplicated(exp$ALIGNED_SEQ))
s.exp3 <- subset(exp, !duplicated(exp$RECONSTRUCTED_SEQ))
s.iw7.del <- subset(iw7.del, !duplicated(iw7.del$RECONSTRUCTED_SEQ))
# for some reason it missed...but it shows up as a second alignment....and i dont know why. It ultimately doesnt matter
#Just remove the duplicates from the merge
s.ip.merge <- subset(ip.merge, !duplicated(ip.merge$RECONSTRUCTED_SEQ))
# YAY it works! Now figure out how to do weighted counts or something
max.ip <- sum(s.ip.merge$percent)
abj.i <- subset(s.ip.merge, REPAIR_TYPE=="ABJ")
abj.ip <- sum(abj.i$percent)
rel.abj <- abj.ip/max.ip*100 # 36.3% of Deletion junctions are ABJ
mhj.i <- subset(s.ip.merge, REPAIR_TYPE=="MHJ")
mhj.ip <- sum(mhj.i$percent)
rel.mhj <- mhj.ip/max.ip*100 # 63.7% of Deletion junctions are MHJ
abj.c <- sum(abj.i[which(abj.i[,30]=="TRUE"),24])
abj.c/abj.ip*100
#[1] 88.17853
mhj.c <- sum(mhj.i[which(mhj.i[,30]=="TRUE"),24])
mhj.c/mhj.ip*100
#[1] 65.61895
