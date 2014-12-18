Signature Validations
-------------------------------------------------------------------------------------------------------
This code corelates GBM and LUAD EGFR mutation/amplification/RPPA data with EGFR signature predictions.

First work with GBM, then LUAD 

library(gmodels)
Load the GBM probabilitity file

binRegResultsGBM <- read.delim("~/Documents/ThesisWork/U01Project/binregResults/EGFR/EGFR_GBM_1Metagene_dwd_200/probabilities.txt", quote="")
head(binRegResultsGBM)

Clean Up the Data: change the sample names from factors to charcters, Remove the contols, and shorten IDS

GBMsamples = binRegResultsGBM$Sample
length(GBMsamples) #206
GBMcharacterSamples= as.character(GBMsamples)
GBMTCGAonly=GBMcharacterSamples[37:187]
length(GBMTCGAonly)
GBMprobabilites=binRegResultsGBM$Probability[37:187]
length(GBMprobabilites)
hist(GBMprobabilites)
GBMshortIDs=substring(GBMTCGAonly, 1, 12)
GBMshortIDs
length(GBMshortIDs)
GBMshortIDs=as.factor(GBMshortIDs)
class(GBMshortIDs)

Make a new data frame with only probs and sampleIDs

GBMprobTable= data.frame(GBMshortIDs,GBMprobabilites)
GBMprobTable

Bring in the mutation status data from Steves File

#Read in the EGFR Class File
EGFRStatus <- read.delim("~/Documents/ThesisWork/U01Project/SignatureValidations/EGFR/PANCAN12_EGFR_MutationStatus.txt", header=FALSE)
EGFRStatus

#merge the EGFR file with the probabilities File
#GBMmergeMutationsAll=merge(GBMprobTable, EGFRStatus, by=1, all.x=TRUE)
GBMmergeMutations=merge(GBMprobTable, EGFRStatus, by=1, all.x=TRUE)
#GBMmergeMutations=merge(GBMprobTable, EGFRStatus, by=1)
colnames(GBMmergeMutations)[3] = "EGFR_Status"
dim(GBMmergeMutations) #151
GBMmergeMutations
summary(GBMmergeMutations) 


# cant do the corelation because they must be numeric. 
cor(GBMmergeMutations$GBMprobabilites,GBMmergeMutations$EGFR_Status)
# box plot the mutation status
plot(GBMmergeMutations$EGFR_Status, GBMmergeMutations$GBMprobabilites, main= "EGFR", ylab = "Probability")
plot(GBMmergeMutations$GBMprobabilites, main="EGFR Probabilities GBM", ylab="Probabilities", col=as.factor(GBMmergeMutations$EGFR_Status))

# break them into mut and not mutated
GBMmutated=subset(GBMmergeMutations,GBMmergeMutations$EGFR_Status=="Mutated" )
GBMmutated
hist(GBMmutated$GBMprobabilites, main= "EGFR Mutated", xlab="Probability")
plot(GBMmutated$GBMprobabilites,main= "EGFR Mutated", ylab="Probability")
GBMnonmutated=subset(GBMmergeMutations,GBMmergeMutations$EGFR_Status=="NotMutated" )
GBMnonmutated
hist(GBMnonmutated$GBMprobabilites, main= "EGFR Not Mutated", xlab="Probability")
plot(GBMnonmutated$GBMprobabilites, main= "EGFR Not Mutated", ylab="Probability")

t.test(GBMmutated$GBMprobabilites, GBMnonmutated$GBMprobabilites)#0.463
wilcox.test(GBMmutated$GBMprobabilites, GBMnonmutated$GBMprobabilites,  conf.int =TRUE) #0.33

#########################################################
# Bring in the amplification/mut data together from cBioPortal. This is the case matrix, which set
# the amp. threshold at over 2. 
#########################################################

CBioPortal_GBM_CaseMatrix_EGFR_AMP_Muts_RNAseq <- read.delim("~/Documents/ThesisWork/U01Project/SignatureValidations/EGFR/GBM/CBioPortal_GBM_CaseMatrix_EGFR_AMP_Muts_RNAseq.txt", header=FALSE)
head(CBioPortal_GBM_CaseMatrix_EGFR_AMP_Muts_RNAseq)
dim(CBioPortal_GBM_CaseMatrix_EGFR_AMP_Muts_RNAseq)
GBM_AMP_Mut_Status=CBioPortal_GBM_CaseMatrix_EGFR_AMP_Muts_RNAseq[3:156,1:2]
head(GBM_AMP_Mut_Status)
dim(GBM_AMP_Mut_Status)
colnames(GBM_AMP_Mut_Status)[2] = "EGFR_Amp_or_Mut"
colnames(GBM_AMP_Mut_Status)[1] = "TCGA_ID"
colnames(GBM_AMP_Mut_Status)

#merge with GBM proabilities/mut file
GBMmergeMutation_AmpsALL=merge(GBMmergeMutations, GBM_AMP_Mut_Status, by=1, all.x=TRUE)
GBMmergeMutation_AmpsALL
summary(GBMmergeMutation_AmpsALL)
table(GBMmergeMutation_AmpsALL$EGFR_Amp_or_Mut) 0= 66, 1= 77
dim(GBMmergeMutation_AmpsALL) #151
GBMmergeMutation_Amps=merge(GBMmergeMutations, GBM_AMP_Mut_Status, by=1)
View(GBMmergeMutation_Amps)
dim(GBMmergeMutation_Amps) #143 removed the NAs so I could do a correlation

cor(GBMmergeMutation_Amps$GBMprobabilites,GBMmergeMutation_Amps$EGFR_Amp_or_Mut) #-0.083

# box plot the mutation$EGFR_Status, GBMmergeMutations$GBMprobabilites, main= "EGFR", ylab = "Probability")
#boxplot(GBMmergeMutation_Amps$GBMprobabilites, GBMmergeMutation_Amps$EGFR_Amp_or_Mut, main= "EGFR")
plot(GBMmergeMutation_Amps$GBMprobabilites, GBMmergeMutation_Amps$EGFR_Amp_or_Mut)

# break them into 1(altered) and 1(not altered)
GBM_altered=subset(GBMmergeMutation_Amps,GBMmergeMutation_Amps$EGFR_Amp_or_Mut=="1" )
GBM_altered
dim(GBM_altered)
summary(GBM_altered)
table(GBM_altered$EGFR_Amp_or_Mut)
hist(GBM_altered$GBMprobabilites, main= "EGFR Altered", xlab="Probability")
plot(GBM_altered$GBMprobabilites,main= "EGFR Altered", ylab="Probability")
GBM_unaltered=subset(GBMmergeMutation_Amps,GBMmergeMutation_Amps$EGFR_Amp_or_Mut=="0" )
GBM_unaltered 
dim(GBM_unaltered)
hist(GBM_unaltered$GBMprobabilites, main= "EGFR Normal", xlab="Probability")
plot(GBM_unaltered$GBMprobabilites,main= "EGFR Normal", ylab="Probability")

boxplot(GBM_unaltered$GBMprobabilites, GBM_altered$GBMprobabilites, main= "Unaltered vs. Altered")
t.test(GBM_altered$GBMprobabilites, GBM_unaltered$GBMprobabilites) #0.319
wilcox.test(GBM_altered$GBMprobabilites, GBM_unaltered$GBMprobabilites) # 0.2896

# Now try with the GISTIC scores, this is going to include the 1s and 2s, but just the 2s
TCGA_GISTIC_GBM_AllTumors <- read.delim("~/Documents/ThesisWork/U01Project/SignatureValidations/EGFR/GBM/GISTIC_GBM_EGFR_RNA-seqOnly.txt", header=FALSE)
head(TCGA_GISTIC_GBM_AllTumors)
dim(TCGA_GISTIC_GBM_AllTumors) #156
# remove the header
GBM_GISTIC=TCGA_GISTIC_GBM_AllTumors[3:156,1:2]
dim(GBM_GISTIC) #154
colnames(GBM_AMP_Mut_Status)[2] = "EGFR_Amp_or_Mut"
colnames(GBM_AMP_Mut_Status)

#merge with GBM proabilities/mut file
GBMmergeMutation_GISTIC_all=merge(GBMmergeMutations, GBM_GISTIC, by=1, all.x=TRUE)
GBMmergeMutation_GISTIC_all
dim(GBMmergeMutation_GISTIC_all) #151
summary(GBMmergeMutation_GISTIC_all) 
# There are alot of 1s. 0=14, 1=58, 2=66
GBMmergeMutation_GISTIC=merge(GBMmergeMutations, GBM_GISTIC, by=1)
GBMmergeMutation_GISTIC
dim(GBMmergeMutation_GISTIC) #143 
colnames(GBMmergeMutation_GISTIC)[4] = "EGFR_Amp_or_Mut"
head(GBMmergeMutation_GISTIC)

#Learning how to make frequency tables
summary(GBMmergeMutation_GISTIC)
CrossTable(GBMmergeMutation_GISTIC$EGFR_Status, GBMmergeMutation_GISTIC$EGFR_Amp_or_Mut)
table(GBMmergeMutation_GISTIC$EGFR_Status)
# gives the frequencies 
table(GBMmergeMutation_GISTIC$EGFR_Amp_or_Mut) #0=14, 1=58, 2=66
#give percents and requires a table
prop.table(table(GBMmergeMutation_GISTIC$EGFR_Amp_or_Mut))
cor(GBMmergeMutation_GISTIC$GBMprobabilites,GBMmergeMutation_GISTIC$EGFR_Amp_or_Mut) 
class(GBMmergeMutation_GISTIC$EGFR_Amp_or_Mut)
*****Change to a number*****
# box plot the mutation status
plot(GBMmergeMutation_GISTIC$EGFR_Status, GBMmergeMutation_GISTIC$GBMprobabilites, main= "EGFR", ylab = "Probability")
plot(GBMmergeMutation_GISTIC$GBMprobabilites, main="EGFR Probabilities GBM", ylab="Probabilities")

# break them into 1(altered) and 1(not altered)
#GBMmergeMutation_GISTIC$EGFR_Amp_or_Mut <- as.numeric(GBMmergeMutation_GISTIC$EGFR_Amp_or_Mut)
GBM_altered=subset(GBMmergeMutation_GISTIC,GBMmergeMutation_GISTIC$EGFR_Amp_or_Mut!="0")
head(GBM_altered)
summary(GBM_altered)
#The summary shows that there are more 1s than there are mutated samples so that is going to be bad and 
# skew our resuts, so the 1s are going to be too much
hist(GBM_altered$GBMprobabilites, main= "EGFR Altered", xlab="Probability")
plot(GBM_altered$GBMprobabilites,main= "EGFR Altered", ylab="Probability")

GBM_unaltered0=subset(GBMmergeMutation_GISTIC,GBMmergeMutation_GISTIC$EGFR_Amp_or_Mut=="0" )
GBM_unaltered0 
dim(GBM_unaltered0)
hist(GBM_unaltered0$GBMprobabilites, main= "EGFR Normal", xlab="Probability")
plot(GBM_unaltered0$GBMprobabilites,main= "EGFR Normal", ylab="Probability")

t.test(GBM_altered$GBMprobabilites, GBM_unaltered$GBMprobabilites) #0.53
wilcox.test(GBM_altered$GBMprobabilites, GBM_unaltered$GBMprobabilites) # 0.84

#compare the 0 to the 2 for GISTIC

GBM_altered2=subset(GBMmergeMutation_GISTIC,GBMmergeMutation_GISTIC$EGFR_Amp_or_Mut=="2")
head(GBM_altered2)
summary(GBM_altered2)
hist(GBM_altered2$GBMprobabilites, main= "EGFR Altered2", xlab="Probability")
plot(GBM_altered2$GBMprobabilites,main= "EGFR Altered2", ylab="Probability")

GBM_unaltered=subset(GBMmergeMutation_GISTIC,GBMmergeMutation_GISTIC$EGFR_Amp_or_Mut=="0" )
GBM_unaltered 

# T-test is better between 0 and 1
t.test(GBM_altered2$GBMprobabilites, GBM_unaltered0$GBMprobabilites) #0.7096
wilcox.test(GBM_altered2$GBMprobabilites, GBM_unaltered0$GBMprobabilites) # 0.6349

#########################################################
# Proteomics2 - now we are using the TCPA(MD Anderson) Data not TCGA
#########################################################

#Read in the proteomics data and sub set EGFR proteins
proteinDataTCPA <- read.csv("~/Documents/ThesisWork/U01Project/SignatureValidations/EGFR/GBM/TCGA-GBM-L3-S42.csv", row.names=1)
head(proteinDataTCPA) #  215 Samples GBM RPPA
EGFRsamples=subset(proteinDataTCPA, select= c( "EGFR", "EGFR_pY1068", "EGFR_pY1173"))
head(EGFRsamples)
summary(GBMmergeMutation_AmpsALL)
class(EGFRsamples)
dim(EGFRsamples) #Still 215

#merge the data with the RPPA data, only keeps the samples we have protein data for
mergeProtein=NULL
class(GBMmergeMutation_Amps)
dim(GBMmergeMutation_Amps) #151 only had RNA-seq data for 151 samples
mergeProteinAll=merge(GBMmergeMutation_Amps, EGFRsamples, by.x="GBMshortIDs", by.y=0, all.x=TRUE)
mergeProteinAll #151 (They only did RPPA on 151 samples)
summary(mergeProteinAll$EGFR_Status)

#remove the NAs
GBMmergeProtein=merge(GBMmergeMutation_Amps, EGFRsamples, by.x="GBMshortIDs", by.y=0)
GBMmergeProtein
dim(GBMmergeProtein)
summary(GBMmergeProtein$EGFR_Status) # 21 mutated and 39 not-mutated

# Correlate the EGFR protein levels with probabilites
--------------------------------------------------------
# basic stats for all EGFR 
quantile(GBMmergeProtein$EGFR) #protein
quantile(GBMmergeProtein$EGFR_pY1068) # 1163
quantile(GBMmergeProtein$EGFR_pY1173) # 1173

# Correlating each EGFR to probabilties
# Results : Not good low neg. correlation
print("EGFR protein Levels")
cor(GBMmergeProtein$GBMprobabilites,GBMmergeProtein$EGFR) #-0.134723
cor(GBMmergeProtein$GBMprobabilites,GBMmergeProtein$EGFR_pY1068) #-0.2006837
cor(GBMmergeProtein$GBMprobabilites,GBMmergeProtein$EGFR_pY1173) #-0.226

#Now plot the EGFR data to see the distrubtion
plot(GBMmergeProtein$GBMprobabilites,GBMmergeProtein$EGFR, main= "EGFR Protein GBM ", ylab= "EGFR Protein Levels", xlab="Probability")
plot(GBMmergeProtein$GBMprobabilites,GBMmergeProtein$EGFR_pY1068, main= "EGFR1068 Phosphorylation GBM ", ylab= "EGFR Phosphorylation Levels", xlab="Probability")
plot(GBMmergeProtein$GBMprobabilites,GBMmergeProtein$EGFR_pY1173, main= "EGFR1173 Phosphorylation GBM ", ylab= "EGFR Phosphorylation Levels", xlab="Probability")

# histogram shows alot more samples with low EGFR protein levels, which we expect
hist(GBMmergeProtein$EGFR,  main= "EGFR Protein GBM ",xlab= "EGFR Protein Levels")
#histogram shows alot more samples with low EGFR phospho levels, which we expect
hist(GBMmergeProtein$EGFR_pY1068,  main= "EGFR1068 Phosphorylation GBM ",xlab= "EGFR Phosphorylation Levels")
hist(GBMmergeProtein$EGFR_pY1173,  main= "EGFR1173 Phosphorylation GBM ",xlab= "EGFR Phosphorylation Levels")


#correlate EGFR protein data to mutation status, Not the useful to use but still makes sense
# Results: All correlate Mutatation status with protein data and are significant
head(GBMmergeProtein)
plot(GBMmergeProtein$EGFR_Status, GBMmergeProtein$EGFR, main= "EGFR protein level GBM", ylab="EGFR Protein Level")
t.test(subset(GBMmergeProtein$EGFR,GBMmergeProtein$EGFR_Status=="Mutated"),subset(mergeProtein$EGFR,GBMmergeProtein$EGFR_Status=="NotMutated"))
plot(GBMmergeProtein$EGFR_Status, GBMmergeProtein$EGFR_pY1068, main= "EGFR phospho level, GBM", ylab="EGFR Ph Level")
t.test(subset(GBMmergeProtein$EGFR_pY1068,GBMmergeProtein$EGFR_Status=="Mutated"),subset(mergeProtein$EGFR_pY1068,GBMmergeProtein$EGFR_Status=="NotMutated"))
plot(GBMmergeProtein$EGFR_Status, GBMmergeProtein$EGFR_pY1173, main= "EGFR phospho level GBM", ylab="EGFR Protein Level")
t.test(subset(GBMmergeProtein$EGFR_pY1173,GBMmergeProtein$EGFR_Status=="Mutated"),subset(GBMmergeProtein[,6],GBMmergeProtein$EGFR_Status=="NotMutated"))

head(GBMmergeProtein)
# Results: All correlate Mutatation/amp status with protein data and are significant

plot(GBMmergeProtein$EGFR_Amp_or_Mut, GBMmergeProtein$EGFR, main= "EGFR protein level GBM", ylab="EGFR Protein Level")
t.test(subset(GBMmergeProtein$EGFR,GBMmergeProtein$EGFR_Amp_or_Mut=="0"),subset(GBMmergeProtein$EGFR,GBMmergeProtein$EGFR_Amp_or_Mut=="1"))

plot(GBMmergeProtein$EGFR_Amp_or_Mut, mergeProtein$EGFR_pY1068, main= "EGFR phospho level 1068 GBM", ylab="EGFR Ph Level")
t.test(subset(GBMmergeProtein$EGFR_pY1068,GBMmergeProtein$EGFR_Status=="0"),subset(GBMmergeProtein$EGFR_pY1068,GBMmergeProtein$EGFR_Status=="1"))

plot(GBMmergeProtein$EGFR_Amp_or_Mut, GBMmergeProtein$EGFR_pY1173, main= "EGFR phospho level 1173 GBM", ylab="EGFR Protein Level")
t.test(subset(GBMmergeProtein$EGFR_pY1173,GBMmergeProtein$EGFR_Status=="1"),subset(GBMmergeProtein[,6],GBMmergeProtein$EGFR_Status=="1"))
hist(GBMmergeProtein$EGFR)
# Try to filter based on the median which equals 0 I belive
EGFR_high= subset(GBMmergeProtein,(GBMmergeProtein$EGFR> 0))
EGFR_high
EGFR_Low=subset(GBMmergeProtein,(GBMmergeProtein$EGFR< 0))
EGFR_Low
dim(EGFR_Low)


boxplot(EGFR_high$GBMprobabilites, EGFR_Low$GBMprobabilites, main= "EGFR GBM", ylab="Probability")
t.test(EGFR_Low$GBMprobabilites, EGFR_high$GBMprobabilites) #0.123
wilcox.test(EGFR_Low$GBMprobabilites, EGFR_high$GBMprobabilites) #0.022



#########################################################
# Now try to correlate probabilites with LUAD
#########################################################

binRegResultsLUAD <- read.delim("~/Documents/ThesisWork/U01Project/binregResults/EGFR/EGFR_LUAD_1Metagene_dwd_200/probabilities.txt", quote="") 
head(binRegResultsLUAD)

# Change samples to characters
LUADsamples = binRegResultsLUAD$Sample
LUADsamples
LUADcharacterSamples= as.character(LUADsamples)
length(LUADcharacterSamples)

#remove the contols
LUADTCGAonly=LUADcharacterSamples[37:531]
length(LUADTCGAonly)
LUADprobabilites=binRegResultsLUAD$Probability[37:531]
length(LUADprobabilites)

#Shorten the TCGA IDs
LUADshortIDs=substring(LUADTCGAonly, 1, 12)
LUADshortIDs
length(LUADshortIDs)
LUADshortIDs=as.factor(LUADshortIDs)
class(LUADshortIDs)

#make a new data frame with only probs and sampleIDs
LUADprobTable= data.frame(LUADshortIDs,LUADprobabilites)
class(LUADprobTable)
LUADprobTable[1,]
#Bimodal distribution
hist(LUADprobTable$LUADprobabilites, main="LUAD probabilities")

#merge the EGFR file with the probabilities File
LUADmergeMutations=merge(LUADprobTable, EGFRStatus, by=1, all.x=TRUE)
LUADmergeMutations
summary(LUADmergeMutations)
LUADmergeMutations=merge(LUADprobTable, EGFRStatus, by=1)
summary(LUADmergeMutations)
LUADmergeMutations
colnames(LUADmergeMutations)[3] = "EGFR_Status"
head(LUADmergeMutations)

# box plot the mutation status
# it looks like the non-mutated samples have higher probabilies :( I think binreg might be too sensitive
#, or that these samples are acutally have the EGFR pathway on but no mutations. 

plot(LUADmergeMutations$EGFR_Status,LUADmergeMutations$LUADprobabilites, main= "EGFR", ylab = "Probability")

plot(LUADmergeMutations$LUADprobabilites, main="EGFR Probabilities GBM", ylab="Probabilities")

# break them into mut and not mutated
# There are more muated samples with probs below 0.5 than over 0.5. 
LUADmutated=subset(LUADmergeMutations,LUADmergeMutations$EGFR_Status=="Mutated" )
LUADmutated # 29
summary(LUADmutated)
hist(LUADmutated$LUADprobabilites, main= "EGFR Mutated", xlab="Probability")
plot(LUADmutated$LUADprobabilites,main= "EGFR Mutated", ylab="Probability")

LUADnonmutated=subset(LUADmergeMutations,LUADmergeMutations$EGFR_Status=="NotMutated" )
LUADnonmutated #166
#bimodal 
hist(LUADnonmutated$LUADprobabilites, main= "EGFR Not Mutated LUAD", xlab="Probability")
# looks like the unmutated have alot of high probabilites. Is this 
plot(LUADnonmutated$LUADprobabilites, main= "EGFR Not Mutated LUAD", ylab="Probability")
t.test(LUADmutated$LUADprobabilites, LUADnonmutated$LUADprobabilites) #0.84
wilcox.test(LUADmutated$LUADprobabilites, LUADnonmutated$LUADprobabilites) #0.876

#########################################################
# Bring in the amplification/mut data together from cBioPortal for LUAD
#########################################################

CBioPortal_LUAD_CaseMatrix_EGFR_AMP_Muts_RNAseq <- read.delim("~/Documents/ThesisWork/U01Project/SignatureValidations/EGFR/LUAD/CBioPortal_LUAD_CaseMatrix_EGFR_AMP_MUTS_RNAseq.txt", header=FALSE )
LUAD_AMP_Mut_Status=CBioPortal_LUAD_CaseMatrix_EGFR_AMP_Muts_RNAseq
head(LUAD_AMP_Mut_Status)
dim(LUAD_AMP_Mut_Status)
colnames(LUAD_AMP_Mut_Status)[2] = "EGFR_Amp_or_Mut"
colnames(LUAD_AMP_Mut_Status)[1] = "TCGA_ID"
head(LUAD_AMP_Mut_Status)

#merge with GBM proabilities/mut file
dim(LUADmergeMutations)
LUADmergeMutation_AmpsALL=merge(LUADmergeMutations, LUAD_AMP_Mut_Status, by=1, all.x=TRUE)
LUADmergeMutation_AmpsALL
summary(LUADmergeMutation_AmpsALL)
dim(LUADmergeMutation_AmpsALL)
LUADmergeMutation_Amps=merge(LUADmergeMutations, LUAD_AMP_Mut_Status, by=1)
  LUADmergeMutation_Amps
dim(LUADmergeMutation_Amps)
summary(LUADmergeMutation_Amps)

cor(LUADmergeMutation_Amps$LUADprobabilites,LUADmergeMutation_Amps$EGFR_Amp_or_Mut) #0.0802

# box plot the mutation status
plot(LUADmergeMutations$EGFR_Status, LUADmergeMutations$GBMprobabilites, main= "EGFR", ylab = "Probability")
plot(LUADmergeMutations$LUADprobabilites, main="EGFR Probabilities LUAD", ylab="Probabilities")

# break them into 1(altered) and 1(not altered)
LUAD_altered=subset(LUADmergeMutation_Amps,LUADmergeMutation_Amps$EGFR_Amp_or_Mut=="1" )
LUAD_altered
dim(LUAD_altered)
hist(LUAD_altered$LUADprobabilites, main= "EGFR Altered", xlab="Probability")
plot(LUAD_altered$LUADprobabilites,main= "EGFR Altered", ylab="Probability")
LUAD_unaltered=subset(LUADmergeMutation_Amps,LUADmergeMutation_Amps$EGFR_Amp_or_Mut=="0" )
LUAD_unaltered 
dim(LUAD_unaltered)
hist(LUAD_unaltered$LUADprobabilites, main= "EGFR Normal", xlab="Probability")
plot(LUAD_unaltered$LUADprobabilites,main= "EGFR Normal", ylab="Probability")

t.test(LUAD_altered$LUADprobabilites, LUAD_unaltered$LUADprobabilites) #0.340

wilcox.test(LUAD_altered$LUADprobabilites, LUAD_unaltered$LUADprobabilites) # .3112

-----------------
LUAD Proteomics
-----------------
#########################################################
# Proteomics2 - now we are using the TCPA(MD Anderson) Data not TCGA
#########################################################

#Read in the LUAD proteomics data and sub set EGFR proteins
LUADproteinDataTCPA <- read.csv("~/Documents/ThesisWork/U01Project/SignatureValidations/EGFR/LUAD/TCGA-LUAD-L3-S51.csv", row.names=1)
View(LUADproteinDataTCPA)
dim(LUADproteinDataTCPA) #  237
LUADEGFRsamples=subset(LUADproteinDataTCPA, select= c( "EGFR", "EGFR_pY1068", "EGFR_pY1173"))
head(LUADEGFRsamples)
dim(LUADEGFRsamples) #147

#merge the data with the RPPA data, only keeps the samples we have protein data for
dim(LUADmergeMutation_Amps) #147 only had RNA-seq data for 147 samples
head(EGFRsamples)
head(LUADmergeMutation_Amps)
rownames(LUADmergeMutation_Amps)= LUADmergeMutation_Amps$LUADshortIDs

LUADmergeProteinAll=merge(LUADmergeMutation_Amps, LUADEGFRsamples, by.x="LUADshortIDs", by.y=0, all.x=TRUE)
summary(LUADmergeProteinAll)
dim(LUADmergeProteinAll)

#remove the NAs
LUADmergeProtein=merge(LUADmergeMutation_Amps, LUADEGFRsamples, by.x="LUADshortIDs", by.y=0)
head(LUADmergeProtein)
dim(LUADmergeProtein) #113
summary(LUADmergeProtein$EGFR_Status) # 19 mutated and 94 not-mutated

  # basic stats for all EGFR 
quantile(LUADmergeProtein$EGFR) #protein
quantile(LUADmergeProtein$EGFR_pY1068) # 1163
quantile(LUADmergeProtein$EGFR_pY1173) # 1173

# Correlating each EGFR to probabilties
# Results : Not good low neg. correlation
print("EGFR protein Levels")
LUADmergeProtein$EGFR
cor(LUADmergeProtein$LUADprobabilites,LUADmergeProtein$EGFR) #-0.03
cor(LUADmergeProtein$LUADprobabilites,LUADmergeProtein$EGFR_pY1068) #0.012
cor(LUADmergeProtein$LUADprobabilites,LUADmergeProtein$EGFR_pY1173) #-0.010

#Now plot the EGFR data to see the distrubtion
plot(LUADmergeProtein$LUADprobabilites,LUADmergeProtein$EGFR, main= "EGFR Protein LUAD ", ylab= "EGFR Protein Levels", xlab="Probability")
plot(LUADmergeProtein$LUADprobabilites,LUADmergeProtein$EGFR_pY1068, main= "EGFR1068 Phosphorylation LUAD ", ylab= "EGFR Phosphorylation Levels", xlab="Probability")
plot(LUADmergeProtein$LUADprobabilites,LUADmergeProtein$EGFR_pY1173, main= "EGFR1173 Phosphorylation LUAD ", ylab= "EGFR Phosphorylation Levels", xlab="Probability")

# histogram shows alot more samples with low EGFR protein levels, which we expect
hist(LUADmergeProtein$EGFR,  main= "EGFR Protein LUAD ",xlab= "EGFR Protein Levels")
#histogram shows alot more samples with low EGFR phospho levels, which we expect
hist(LUADmergeProtein$EGFR_pY1068,  main= "EGFR1068 Phosphorylation LUAD ",xlab= "EGFR Phosphorylation Levels")
hist(LUADmergeProtein$EGFR_pY1173,  main= "EGFR1173 Phosphorylation LUAD ",xlab= "EGFR Phosphorylation Levels")


#correlate EGFR protein data to mutation status, Not the useful to use but still makes sense
# Results: All correlate Mutatation status with protein data and are significant
head(LUADmergeProtein)
plot(LUADmergeProtein$EGFR_Status, LUADmergeProtein$EGFR, main= "EGFR protein level LUAD", ylab="EGFR Protein Level")
t.test(subset(LUADmergeProtein$EGFR,LUADmergeProtein$EGFR_Status=="Mutated"),subset(LUADmergeProtein$EGFR,LUADmergeProtein$EGFR_Status=="NotMutated"))
plot(LUADmergeProtein$EGFR_Status, LUADmergeProtein$EGFR_pY1068, main= "EGFR phospho level, LUAD", ylab="EGFR Ph Level")
t.test(subset(LUADmergeProtein$EGFR_pY1068,LUADmergeProtein$EGFR_Status=="Mutated"),subset(LUADmergeProtein$EGFR_pY1068,LUADmergeProtein$EGFR_Status=="NotMutated"))
plot(LUADmergeProtein$EGFR_Status, LUADmergeProtein$EGFR_pY1173, main= "EGFR phospho level LUAD", ylab="EGFR Protein Level")
t.test(subset(LUADmergeProtein$EGFR_pY1173,LUADmergeProtein$EGFR_Status=="Mutated"),subset(LUADmergeProtein[,6],LUADmergeProtein$EGFR_Status=="NotMutated"))

head(GBMmergeProtein)
# Results: All correlate Mutatation/amp status with protein data and are significant

plot(LUADmergeProtein$EGFR_Amp_or_Mut, LUADmergeProtein$EGFR, main= "EGFR protein level LUAD", ylab="EGFR Protein Level")
t.test(subset(LUADmergeProtein$EGFR,LUADmergeProtein$EGFR_Amp_or_Mut=="0"),subset(LUADmergeProtein$EGFR,LUADmergeProtein$EGFR_Amp_or_Mut=="1"))

plot(LUADmergeProtein$EGFR_Amp_or_Mut, LUADmergeProtein$EGFR_pY1068, main= "EGFR phospho level 1068 LUAD", ylab="EGFR Ph Level")
t.test(subset(GBMmergeProtein$EGFR_pY1068,GBMmergeProtein$EGFR_Status=="0"),subset(GBMmergeProtein$EGFR_pY1068,GBMmergeProtein$EGFR_Status=="1"))

plot(LUADmergeProtein$EGFR_Amp_or_Mut, LUADmergeProtein$EGFR_pY1173, main= "EGFR phospho level 1173 GBM", ylab="EGFR Protein Level")
t.test(subset(LUADmergeProtein$EGFR_pY1173,LUADmergeProtein$EGFR_Status=="1"),subset(LUADmergeProtein[,6],LUADmergeProtein$EGFR_Status=="1"))

# Try to filter based on the median which equals 0 I belive
LUAD_EGFR_high= subset(LUADmergeProtein,(LUADmergeProtein$EGFR>0))
LUAD_EGFR_high
dim(LUAD_EGFR_high)
LUAD_EGFR_Low=subset(LUADmergeProtein,(LUADmergeProtein$EGFR<0))
LUAD_EGFR_Low
dim(LUAD_EGFR_Low)

boxplot(LUAD_EGFR_high$LUADprobabilites, LUAD_EGFR_Low$LUADprobabilites, main= "EGFR LUAD", ylab="Probability")
t.test(LUAD_EGFR_Low$LUADprobabilites, LUAD_EGFR_high$LUADprobabilites) #0.1132
wilcox.test(LUAD_EGFR_Low$LUADprobabilites, LUAD_EGFR_high$LUADprobabilites) #0.1054








#########################################################
# Proteomics- This was Using the TCGA data, not doing to 
# use this anymore
#########################################################

#Read in the proteomics data
#proteinData <- t(read.delim('/Users/Shasta/Desktop/ThesisWork/mdanderson.org_GBM_MDA_RPPA_Core.txt', quote="", check.names=FALSE, header = 1, row.names= 1))
#class(proteinData)

#subet EGFR data 
#EGFRsamples=subset(proteinData, select= c("EGFR-R-C", "EGFR_pY1068-R-V", "EGFR_pY1173-R-C", "EGFR_pY992-R-V"))
#EGFRsamples

#dim(EGFRsamples)

#remove contols and shorten IDs
#EGFRsamples2=EGFRsamples[grepl("TCGA*", rownames(EGFRsamples)), ]
#dim(EGFRsamples)
#EGFRsamples2
#proteinShortIDs=substring(rownames(EGFRsamples2), 1, 12)
#proteinShortIDs=as.factor(proteinShortIDs)
#row.names(EGFRsamples2)= proteinShortIDs

#merge the data with the RPPA data, only keeps the samples we have protein data for
#mergeProtein=NULL
#mergeProtein=merge(mergeMutations, EGFRsamples2, by.x="shortIDs", by.y=0)
#print( "21 Mutated, 45 Not")

# Correlate the EGFR protein levels with probabilites
#########################################################

# Results : Correlations not good and inverse

# basic stats for all EGFR 
#quantile(mergeProtein[,4]) #protein
#quantile(mergeProtein[,5]) # 1163
#quantile(mergeProtein[,6]) # 1173
#quantile(mergeProtein[,7]) #922

# Correlating each EGFR to probabilties
#print("EGFR protein Levels")
#EGFR_cor= cor(mergeProtein$ProbCol,mergeProtein[,4])
# #EGFR_cor #-0.22
# #print("EGFR-Ty1068(activating) ")
# EGFR1068_cor=cor(mergeProtein$ProbCol,mergeProtein[,5])
# EGFR1068_cor #-0.20
# print("EGFR-Ty1173")
# EGFR1173_cor=cor(mergeProtein$ProbCol,mergeProtein[,6])
# EGFR1173_cor #-0.22
# print("EGFR-Ty992")
# EGFR992_cor=cor(mergeProtein$ProbCol,mergeProtein[,7])
# EGFR992_cor #-0.24
# # None are correlating, looks like this is due to to the outliers? 
# 
# #Now plot the EGFR data to see the distrubtion
# plot(mergeProtein$ProbCol,mergeProtein[,4], main= "EGFR Protein GBM ", ylab= "EGFR Protein Levels", xlab="Probability")
# plot(mergeProtein$ProbCol,mergeProtein[,5], main= "EGFR1068 Phosphorylation GBM ", ylab= "EGFR Phosphorylation Levels", xlab="Probability")
# plot(mergeProtein$ProbCol,mergeProtein[,6], main= "EGFR1173 Phosphorylation GBM ", ylab= "EGFR Phosphorylation Levels", xlab="Probability")
# plot(mergeProtein$ProbCol,mergeProtein[,7], main= "EGFR992 Phosphorylation GBM ", ylab= "EGFR Phosphorylation Levels", xlab="Probability")
# 
# # histogram shows alot more samples with low EGFR protein levels, which we expect
# hist(mergeProtein[,4],  main= "EGFR Protein GBM ",xlab= "EGFR Protein Levels")
# #histogram shows alot more samples with low EGFR phospho levels, which we expect
# hist(mergeProtein[,5],  main= "EGFR1068 Phosphorylation GBM ",xlab= "EGFR Phosphorylation Levels")
# hist(mergeProtein[,6],  main= "EGFR1173 Phosphorylation GBM ",xlab= "EGFR Phosphorylation Levels")
# hist(mergeProtein[,7],  main= "EGFR992 Phosphorylation GBM ",xlab= "EGFR Phosphorylation Levels")
# 
# #correlate EGFR protein data to mutation status
# # Results: protein and phospho levels are higher in mutateds samples, significant
# plot(mergeProtein$EGFR_Status, mergeProtein[,4], main= "EGFR protein level GBM", ylab="EGFR Protein Level")
# t.test(subset(mergeProtein[,4],mergeProtein$EGFR_Status=="Mutated"),subset(mergeProtein[,4],mergeProtein$EGFR_Status=="NotMutated"))
# 
# plot(mergeProtein$EGFR_Status, mergeProtein[,5], main= "EGFR phospho level, activating GBM", ylab="EGFR Ph Level")
# t.test(subset(mergeProtein[,5],mergeProtein$EGFR_Status=="Mutated"),subset(mergeProtein[,5],mergeProtein$EGFR_Status=="NotMutated"))
# 
# plot(mergeProtein$EGFR_Status, mergeProtein[,6], main= "EGFR phospho level repressing GBM", ylab="EGFR Protein Level")
# t.test(subset(mergeProtein[,6],mergeProtein$EGFR_Status=="Mutated"),subset(mergeProtein[,6],mergeProtein$EGFR_Status=="NotMutated"))
# 
# plot(mergeProtein$EGFR_Status, mergeProtein[,7], main= "EGFR phospho level GBM", ylab="EGFR Phospho Level")
# t.test(subset(mergeProtein[,7],mergeProtein$EGFR_Status=="Mutated"),subset(mergeProtein[,7],mergeProtein$EGFR_Status=="NotMutated"))
# 
# # Try and correlate with the EGFR protein data
# EGFR_high= subset(mergeProtein,(mergeProtein[,4]>-0.5))
# EGFR_high
# dim(EGFR_high)
# EGFR_Low=subset(mergeProtein,(mergeProtein[,4]<=-0.5))
# EGFR_Low
# dim(EGFR_Low)
# t.test(EGFR_Low$ProbCol, EGFR_high$ProbCol) #0.1095
# EGFR_cor_sep= cor(EGFR_high$ProbCol,EGFR_high[,4]) # -0.219
# EGFR_cor_sep
# 
# EGFR_high= subset(mergeProtein,(mergeProtein[,4]>1))
# EGFR_high
# dim(EGFR_high)
# EGFR_Low=subset(mergeProtein,(mergeProtein[,4]<=1))
# EGFR_Low
# dim(EGFR_Low)
# t.test(EGFR_Low$ProbCol, EGFR_high$ProbCol) #0.11
# EGFR_cor_sep= cor(EGFR_high$ProbCol,EGFR_high[,4]) # -0.15
# EGFR_cor_sep
# 
# # seperate based on -0.5 
# EGFR_high= subset(mergeProtein,(mergeProtein[,5]>-0.5))
# EGFR_high
# dim(EGFR_high)
# EGFR_Low=subset(mergeProtein,(mergeProtein[,5]<=-0.5))
# EGFR_Low
# dim(EGFR_Low)
# t.test(EGFR_Low$ProbCol, EGFR_high$ProbCol) #0.02
# EGFR_cor_sep= cor(EGFR_high$ProbCol,EGFR_high[,5])
# EGFR_cor_sep #0.09
# 
# # try with median
# EGFR_high= subset(mergeProtein,(mergeProtein[,5]>-1.7))
# EGFR_high
# dim(EGFR_high)
# EGFR_Low=subset(mergeProtein,(mergeProtein[,5]<=-1.7))
# EGFR_Low
# dim(EGFR_Low)
# t.test(EGFR_Low$ProbCol, EGFR_high$ProbCol) #0.2
# EGFR_cor_sep= cor(EGFR_high$ProbCol,EGFR_high[,5])
# EGFR_cor_sep #-0.23
# 
# #try with the 25% and 75 %
# EGFR_high= subset(mergeProtein,(mergeProtein[,5]<-2.6))
# EGFR_high
# dim(EGFR_high)
# EGFR_Low=subset(mergeProtein,(mergeProtein[,5]<=-1.7))
# EGFR_Low
# dim(EGFR_Low)
# t.test(EGFR_Low$ProbCol, EGFR_high$ProbCol) #0.2
# EGFR_cor_sep= cor(EGFR_high$ProbCol,EGFR_high[,5])
# EGFR_cor_sep #-0.23



