args=(commandArgs(TRUE))
outDir = commandArgs()[7]
library(data.table)

#outDir = "/data2/Moom_datasets"
setwd("/data4/MoomFiles/moom_rahman/PANCAN24")
print(paste("Reading files from",getwd(), sep=" "))

pancan24<-data.frame(fread("06_01_15_TCGA_24.tumor_Rsubread_TPM.txt"),row.names=T,check.names=F)
samples<-read.table("06_01_15_TCGA_24_CancerType_Samples.txt",row.names=1)
Levels=levels(samples$V2)
print(paste("Saving outputs in ",outDir, sep=' '))
for(i in 1:length(Levels)){
    subtype<-subset(pancan24,select=c(samples[colnames(pancan24),1]==Levels[i]))
    outFile<-paste(outDir,paste("PANCAN24",Levels[i],ncol(subtype),"TPMlog2.txt",sep='_'),sep="/")
    print(dim(subtype))
    write.table(log2(subtype+1),outFile,col.names=NA, sep='\t', quote=F)
    }

