short_to_long_TCGA_id=function(longnames=NULL,shortnames=NULL){
  counter=0
  for (j in 1:length(shortnames)){
    if(!is.na(pmatch(shortnames[j],longnames))){
      shortnames[j]<-longnames[pmatch(shortnames[j],longnames, duplicates.ok=F)]  
      counter=counter+1
    }
  }
  print(paste(counter,"names have been changed",sep= " "))
  return(shortnames)
}

tcga_brca_metastasis<-read.table("~/Desktop/TCGA_brca_metastasis_yes.txt", head=T, sep='\t',row.names=1,stringsAsFactors = F)
dim(tcga_brca_metastasis)
rownames(tcga_brca_metastasis)<-gsub("-",".",rownames(tcga_brca_metastasis))
head(tcga_brca_metastasis)
rownames(tcga_brca_metastasis)<-short_to_long_TCGA_id(shortnames = rownames(tcga_brca_metastasis),longnames = rownames(preds))
tcga_brca_no_metastasis<-read.table("~/Desktop/TCGA_brca_metastasis_no.txt", head=T, sep='\t',row.names=1)
dim(tcga_brca_no_metastasis)
rownames(tcga_brca_no_metastasis)<-gsub("-",".",rownames(tcga_brca_no_metastasis))
head(tcga_brca_no_metastasis)
rownames(tcga_brca_no_metastasis)<-short_to_long_TCGA_id(shortnames = rownames(tcga_brca_no_metastasis),longnames = rownames(preds))

preds<-read.csv("~/Dropbox/bild_signatures//tcga_15_mar_all/akt_bad_igf1r_erk/adap_adap_multi/pathway_activity_testset.csv", row.names=1, header = 1)
head(preds)


pam50<-read.table("~/Dropbox/Datasets/tcga_breast_pam50.txt",sep='\t', stringsAsFactors = T,header=T, row.names=1)
partial_sample_names<-rownames(pam50)
sample_names<-rownames(preds)
counter=0
for (j in 1:length(partial_sample_names)){
  if(!is.na(pmatch(partial_sample_names[j],sample_names))){
    partial_sample_names[j]<-sample_names[pmatch(partial_sample_names[j],sample_names, duplicates.ok=F)]  
    counter=counter+1
  }
}

counter
rownames(pam50)<-short_to_long_TCGA_id(longnames = rownames(preds),shortnames = rownames(pam50))
head(pam50)

my_palette <- colorRampPalette(c("darkblue","aliceblue","brown4"))(n = 299)
col_breaks = c(seq(0,0.2,length=100), seq(0.2,0.4,length=100), seq(0.4,1,length=100)) 

met<-merge_drop(preds,tcga_brca_metastasis, by=0)
no_met<-merge_drop(preds,tcga_brca_no_metastasis, by=0)
