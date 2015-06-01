## running ASSIGN for raf sing
trainingLabelr<-list(control=list(raf=1:12),raf=13:18)

sub_dir<-paste(basedir,"raf_25_gene_list",sep='/')
dir.create( sub_dir)
assign_easy_multi(trainingData = cbind(c_gfp,c_raf),test=c_test,trainingLabel1 = trainingLabelr,g=25,out_dir_base = sub_dir,single = 1)

sub_dir<-paste(basedir,"raf_50_gene_list",sep='/')
dir.create( sub_dir)
assign_easy_multi(trainingData = cbind(c_gfp,c_raf),test=c_test,trainingLabel1 = trainingLabelr,g=50,out_dir_base = sub_dir,single = 1)


## accumulating the ASSIGN prediction files in this block
setwd("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/")
filenames<-system("ls *gene_list/*/pathway_activity_testset*", intern=TRUE)
filenames

for(i in 1:length(filenames))
{
  f<-read.csv(filenames[i], header=1,row.names=1) ###reading in the filess one at a time
  colnames(f)<-paste(filenames[i],colnames(f),sep='/')
  if(i!=1){
    print(i)
    data_icbp<-cbind(data_icbp,f)
  }
  else{
    data_icbp<-f
  }
}
#head(data_icbp)
#dim(data_icbp)

colnames(data_icbp)<-gsub(pattern = "/pathway_activity_testset.csv",replacement = "",x = colnames(data_icbp))
head(data_icbp)
rownames(data_icbp)[1:7]<-c("184A1","184B5","21MT1","21MT2","21NT","21PT","600MPE")
setwd("~/Dropbox/bild_signatures//Datasets")
drugs<-read.delim("ICBP_drugs.txt", header=1, sep='\t',row.names=1)
icbp_drug<-merge_drop(data_icbp,drugs)
colnames(icbp_drug)
cor_mat=p_mat=matrix(0,length(filenames),90)
rownames(cor_mat)=rownames(p_mat)=colnames(icbp_drug)[1:length(filenames)]
colnames(cor_mat)=colnames(p_mat)=colnames(icbp_drug)[(length(filenames)+11):ncol(icbp_drug)]

for(i in 1:length(filenames)){
  for(j in 1:90){
    temp=cor.test(icbp_drug[,i],icbp_drug[,(j+length(filenames)+10)],use="pairwise",method="spearman")
    print(j)
    print(temp)
    cor_mat[i,j]=temp$estimate
    p_mat[i,j]=temp$p.value
  }
}
write.table(cor_mat,"single_cor_drug_mat_4_21.txt",col.names = NA,quote=F,sep='\t')
colnames(p_mat)=paste(colnames(p_mat),"p_value",sep="_")
write.table(p_mat,"single_p_drug_mat_4_21.txt",col.names = NA,quote=F,sep='\t')
cor_p_mat<-cbind(cor_mat,p_mat)
order(colnames(cor_p_mat))
cor_p_mat<-cor_p_mat[,order(colnames(cor_p_mat))]
write.table(cor_p_mat,"single_cor_p_mat_4_21.txt",col.names = NA,quote=F,sep='\t')
