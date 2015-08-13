pancan20_tpm_log<-fread("/data/TCGA/PANCAN20/PANCAN20.IlluminaHiSeq_RNASeqV2.tumor_Rsubread_TPMlog_10_9_withclass.txt",header=T,sep='\t')
cancer_type<-function(type=NULL){
    class<-pancan20_tpm_log[pancan20_tpm_log$LGG==type,];
    rownames(class)<-class$V1;
    class_f<-t(subset(class,select=2:(ncol(class)-1)));
    colnames(class_f)<-rownames(class);
    write.table(class_f,paste("/data2/Moom_datasets/TCGA_PANCAN20_Rsubread_",type,"_TPMlog_10_9.txt",sep=""),col.names=NA,sep='\t',quote=F);
    print(head(class_f));
    print(dim(class_f));
    print(type)}
pub<-read.table("/data/TCGA/PANCAN20/TCGA_CancerTypes_Publishable_Sept2014.txt",header=F)
for (i in 1:nrow(pub)){cancer_type(pub[i,1])}

