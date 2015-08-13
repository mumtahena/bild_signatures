library(devtools)
install_github("wevanjohnson/ASSIGN",ref="adapt_gene_only_ver2")
library(ASSIGN)
library(foreach)
library(sva)
library(data.table)
library(doParallel)
library(BatchQC)
cl<-makeCluster(2)
registerDoParallel(cl)


source("~/Dropbox/bild_signatures/bild_signatures/Rmarkdowns_scripts/Key_ASSIGN_functions.Rmd")

setwd("~/Dropbox/bild_signatures/Datasets/")
expr<-as.matrix(read.table("/Users/mumtahenarahman/Dropbox/Datasets/GFP18_AKT_BAD_HER2_IGF1R_RAF_ERK.tpmlog",sep='\t',row.names=1,header=1))
control<-subset(expr, select=GFP.1:GFP.12)
her2<-subset(expr, select=HER2.1:HER2.6)
akt<-subset(expr,select=AKT.1:AKT.6)
bad<-subset(expr,select=BAD.1:BAD.6)
igf1r<-subset(expr,select=IGF1R.1:IGF1R.6)
raf<-subset(expr,select=RAF.1:RAF.6)
erk<-subset(expr,select=ERK.1:ERK.6)
expr_all<-cbind(control,akt,bad,her2,igf1r,raf,erk)
dim(expr_all)
expr_all_f <-expr_all[apply(expr_all[,1:47]==0,1,mean) < 0.85,]
control_egfr_l<-read.table("~/Dropbox/Datasets/18_GFP_EGFR_TPMlog2.txt", sep='\t', header=1, row.names=1)
gfp_egfr_multi_f <- merge_drop(control_egfr_l,expr_all_f)
dim(gfp_egfr_multi_f)
gfp_kras<-read.table("~/Dropbox/Datasets/36_GFP_KRAS_TPMlog2.txt", sep='\t', header=1, row.names=1)
head(gfp_kras)
gfp_egfr_kras_multi_f<-merge_drop(gfp_egfr_multi_f,gfp_kras)
dim(gfp_egfr_kras_multi_f)


icbp<-read.table("~/Dropbox/bild_signatures/Datasets/icbp_Rsubread_tpmlog.txt", sep='\t', header=1, row.names=1)
dim(icbp)
expr_all_icbp_f<-merge_drop(gfp_egfr_kras_multi_f,icbp)
dim(expr_all_icbp_f)
colnames(expr_all_icbp_f)
sub<-c(6,6,12,6,6,5,6,6,6,9,9,9,9,ncol(icbp))
pcaplot(expr_all_icbp_f,sub)
bat<-as.matrix(cbind(colnames(expr_all_icbp_f),c(rep(1,ncol(control_egfr_l)),rep(2,ncol(expr_all)),rep(3,ncol(gfp_kras)),rep(4,ncol(icbp)))))
#modbat<-model.matrix(~as.factor(bat[,2]),data=expr_all_icbp_f)
#batchQC(expr_all_icbp_f, batch=bat1[,2], mod=modbat,report_file="batchqc_report.html", report_dir=".")
combat_expr<-ComBat(dat=expr_all_icbp_f, batch=bat[,2], mod=NULL, numCovs=NULL)
pcaplot(combat_expr,sub)




###########KRAS on TCGA LUAD#########
library(data.table)
gfp_kras<-read.table("~/Dropbox/Datasets/36_GFP_KRAS_TPMlog2.txt", sep='\t', header=1, row.names=1)
head(gfp_kras)
tcga_luad<-data.frame(fread("~/Dropbox/Datasets/PANCAN24_LUAD_541_TPMlog2.txt"),check.names=F,row.names=1)
sub=c(9,9,9,9,541)
expr_f <-gfp_kras[apply(gfp_kras[,1:36]==0,1,mean) < 0.85,]
expr<-merge_drop(expr_f,tcga_luad)
pcaplot(expr,sub)
bat<-as.matrix(cbind(colnames(expr),c(rep(1,ncol(gfp_kras)),rep(2,ncol(tcga_luad)))))
combat_expr<-ComBat(dat=expr, batch=bat[,2], mod=NULL, numCovs=NULL)
pcaplot(combat_expr,sub)

# c_gfp<-subset(combat_expr, select=GFP.1:GFP.12)
# c_akt<-subset(combat_expr, select=AKT.1:AKT.6)
# c_bad<-subset(combat_expr, select=BAD.1:BAD.6)
# c_her2<-subset(combat_expr, select=HER2.1:HER2.6)
# c_igf1r<-subset(combat_expr, select=IGF1R.1:IGF1R.6)
# c_raf<-subset(combat_expr, select=RAF.1:RAF.6)
# c_erk<-subset(combat_expr, select=ERK.1:ERK.6)
# train_egfr<-combat_expr[,1:12]
c_kras_gfp<-subset(combat_expr,select=GFP.31:GFP.39)
c_kraswt<-subset(combat_expr,select=KRASWT.1:KRASWT.9)
c_krasqh<-subset(combat_expr,select=KRASQH.1:KRASQH.9)
c_krasgv<-subset(combat_expr,select=KRASGV.1:KRASGV.9)
c_test<-combat_expr[,(ncol(gfp_kras)+1):ncol(combat_expr)]
colnames(c_test)

#######getting the genelist#######
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/akt_75_gene_list/adapB_single/output.rda")
akt_75_genelist<-output.data$processed.data$diffGeneList  
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/single/akt_150_gene_list/adapB_single/output.rda")
akt_150_genelist<-output.data$processed.data$diffGeneList 
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/bad_200_gene_list/adapB_single/output.rda")
bad_200_genelist<-output.data$processed.data$diffGeneList  
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/igf1r_75_gene_list/adapB_single/output.rda")
igf1r_75_genelist<-output.data$processed.data$diffGeneList  
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/erk_250_gene_list/adapB_single/output.rda")
erk_250_genelist<-output.data$processed.data$diffGeneList  
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/single/erk_400_gene_list/adapB_single/output.rda")
erk_400_genelist<-output.data$processed.data$diffGeneList  
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/her2_15_gene_list/adapB_single/output.rda")
her2_15_genelist<-output.data$processed.data$diffGeneList  
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/egfr_25_gene_list/adapB_single/output.rda")
egfr_25_genelist<-output.data$processed.data$diffGeneList  
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/raf_100_gene_list/adapB_single/output.rda")
raf_100_genelist<-output.data$processed.data$diffGeneList

load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/krasqh_100_gene_list/adapB_single/output.rda")
krasqh_100_genelist<-output.data$processed.data$diffGeneList
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/krasqv_100_gene_list/adapB_single/output.rda")
krasgv_100_genelist<-output.data$processed.data$diffGeneList
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/kraswt_100_gene_list/adapB_single/output.rda")
kraswt_100_genelist<-output.data$processed.data$diffGeneList

load("~/Dropbox/bild_signatures/kras/krasqh_200_gene_list/adapB_single/output.rda")
krasqh_200_genelist<-output.data$processed.data$diffGeneList
load("~/Dropbox/bild_signatures/kras/krasgv_200_gene_list/adapB_single/output.rda")
krasgv_200_genelist<-output.data$processed.data$diffGeneList
load("~/Dropbox/bild_signatures/kras/kraswt_200_gene_list/adapB_single/output.rda")
kraswt_200_genelist<-output.data$processed.data$diffGeneList

load("~/Dropbox/bild_signatures/kras/krasqh_300_gene_list/adapB_single/output.rda")
krasqh_300_genelist<-output.data$processed.data$diffGeneList
load("~/Dropbox/bild_signatures/kras/krasgv_300_gene_list/adapB_single/output.rda")
krasgv_300_genelist<-output.data$processed.data$diffGeneList
load("~/Dropbox/bild_signatures/kras/kraswt_300_gene_list/adapB_single/output.rda")
kraswt_300_genelist<-output.data$processed.data$diffGeneList




#single_pathway_best<-data_icbp[,c("akt_75_gene_list/adap_adap_single/akt","her2_15_gene_list/adap_adap_single/her2","igf1r_75_gene_list/adap_adap_single/igf1r","bad_200_gene_list/adap_adap_single/bad","erk_250_gene_list/adap_adap_single/erk","raf_100_gene_list/adap_adap_single/raf","egfr_25_gene_list/adap_adap_single/egfr")]
#best_pathway_predictions<-pred_drug[,c("akt_75_gene_list/adap_adap_single/akt","her2_25_gene_list/adap_adap_single/her2","igf1r_75_gene_list/adap_adap_single/igf1r","bad_200_gene_list/adap_adap_single/bad","erk_250_gene_list/adap_adap_single/erk","raf_100_gene_list/adap_adap_single/raf","egfr_25_gene_list/adap_adap_single/egfr")]
#single_pathway_best<-c("akt_75_gene_list/adap_adap_single/akt","bad_200_gene_list/adap_adap_single/bad","her2_15_gene_list/adap_adap_single/her2","igf1r_75_gene_list/adap_adap_single/igf1r","erk_200_gene_list/adap_adap_single/erk","raf_100_gene_list/adap_adap_single/raf","egfr_25_gene_list/adap_adap_single/egfr")
#multi_pathway_best<-c("akt_bad_her2_raf/adap_adap_multi/akt","akt_bad_her2_erk/adap_adap_multi/bad","her2_igf1r/adap_adap_multi/her2","igf1r_raf_her2_erk/adap_adap_multi/igf1r","erk_bad/adap_adap_multi/erk","akt_bad_igf1r_raf/adapB_multi/raf")

#####validation in TCGA LUAD 541 samples based on KRAS mutation status. Hypothesis is, mutant sample will have higher kras activity than wild type kras samples;
###KRASWT signature was used as control and expected to see no difference.
basedir="~/Dropbox/bild_signatures/kras_tcga_luad"
dir.create(basedir)
testSig <- function(sigProtein, numGenes=NA, geneList =NULL, trainingData, testData, trainingLabels)
{
  names(sigProtein)<-sigProtein
  trainingLabel<-list(control=list(sigProtein=1:trainingLabels[1]),sigProtein=(trainingLabels[1]+1):(trainingLabels[1]+trainingLabels[2]))
  names(trainingLabel$control)=names(trainingLabel)[2]=names(sigProtein)
  if(is.na(numGenes)){
    sub_dir<-paste(basedir,paste(sigProtein,"gene_list", sep="_"),sep='/')
  }
  else{
    sub_dir<-paste(basedir,paste(sigProtein,numGenes,"gene_list", sep="_"),sep='/')
  }
  dir.create(sub_dir)
  assign_easy_multi(trainingData = trainingData,test=testData,trainingLabel1 = trainingLabel,g=numGenes,geneList = geneList,out_dir_base = sub_dir,single = 1)
}



# You could execute it in parallel like this:
library(doParallel)
library(ASSIGN)
registerDoParallel(cores=2)
trainingLabel<-list(control=list(kraswt=1:9),kraswt=(10:18))
dir.create(paste(basedir,paste("kraswt",100,"gene_list", sep="_"),sep='/'))
trainingLabel<-list(control=list(krasgv=1:9),krasgv=(10:18))
assign_easy_multi(trainingData = cbind(c_kras_gfp,c_kraswt),test=c_test,trainingLabel1 = trainingLabel,geneList = kraswt_100_genelist,out_dir_base = paste(basedir,paste("kraswt",100,"gene_list", sep="_"),sep='/'),single = 1)
trainingLabel<-list(control=list(krasgv=1:9),krasgv=(10:18))
dir.create(paste(basedir,paste("krasgv",100,"gene_list", sep="_"),sep='/'))
trainingLabel<-list(control=list(krasgv=1:9),krasgv=(10:18))
assign_easy_multi(trainingData = cbind(c_kras_gfp,c_krasgv),test=c_test,trainingLabel1 = trainingLabel,geneList = krasgv_100_genelist,out_dir_base = paste(basedir,paste("krasgv",100,"gene_list", sep="_"),sep='/'),single = 1)
trainingLabel<-list(control=list(krasqh=1:9),krasqh=(10:18))
dir.create(paste(basedir,paste("krasqh",100,"gene_list", sep="_"),sep='/'))
assign_easy_multi(trainingData = cbind(c_kras_gfp,c_krasqh),test=c_test,trainingLabel1 = trainingLabel,geneList = krasqh_100_genelist,out_dir_base = paste(basedir,paste("krasqh",100,"gene_list", sep="_"),sep='/'),single = 1)

trainingLabel<-list(control=list(kraswt=1:9),kraswt=(10:18))
dir.create(paste(basedir,paste("kraswt",200,"gene_list", sep="_"),sep='/'))
trainingLabel<-list(control=list(krasgv=1:9),krasgv=(10:18))
assign_easy_multi(trainingData = cbind(c_kras_gfp,c_kraswt),test=c_test,trainingLabel1 = trainingLabel,geneList = kraswt_200_genelist,out_dir_base = paste(basedir,paste("kraswt",200,"gene_list", sep="_"),sep='/'),single = 1)
trainingLabel<-list(control=list(krasgv=1:9),krasgv=(10:18))
dir.create(paste(basedir,paste("krasgv",200,"gene_list", sep="_"),sep='/'))
trainingLabel<-list(control=list(krasgv=1:9),krasgv=(10:18))
assign_easy_multi(trainingData = cbind(c_kras_gfp,c_krasgv),test=c_test,trainingLabel1 = trainingLabel,geneList = krasgv_200_genelist,out_dir_base = paste(basedir,paste("krasgv",200,"gene_list", sep="_"),sep='/'),single = 1)
trainingLabel<-list(control=list(krasqh=1:9),krasqh=(10:18))
dir.create(paste(basedir,paste("krasqh",200,"gene_list", sep="_"),sep='/'))
assign_easy_multi(trainingData = cbind(c_kras_gfp,c_krasqh),test=c_test,trainingLabel1 = trainingLabel,geneList = krasqh_200_genelist,out_dir_base = paste(basedir,paste("krasqh",200,"gene_list", sep="_"),sep='/'),single = 1)


trainingLabel<-list(control=list(kraswt=1:9),kraswt=(10:18))
dir.create(paste(basedir,paste("kraswt",300,"gene_list", sep="_"),sep='/'))
trainingLabel<-list(control=list(krasgv=1:9),krasgv=(10:18))
assign_easy_multi(trainingData = cbind(c_kras_gfp,c_kraswt),test=c_test,trainingLabel1 = trainingLabel,geneList = kraswt_300_genelist,out_dir_base = paste(basedir,paste("kraswt",300,"gene_list", sep="_"),sep='/'),single = 1)
trainingLabel<-list(control=list(krasgv=1:9),krasgv=(10:18))
dir.create(paste(basedir,paste("krasgv",300,"gene_list", sep="_"),sep='/'))
trainingLabel<-list(control=list(krasgv=1:9),krasgv=(10:18))
assign_easy_multi(trainingData = cbind(c_kras_gfp,c_krasgv),test=c_test,trainingLabel1 = trainingLabel,geneList = krasgv_300_genelist,out_dir_base = paste(basedir,paste("krasgv",300,"gene_list", sep="_"),sep='/'),single = 1)
trainingLabel<-list(control=list(krasqh=1:9),krasqh=(10:18))
dir.create(paste(basedir,paste("krasqh",300,"gene_list", sep="_"),sep='/'))
assign_easy_multi(trainingData = cbind(c_kras_gfp,c_krasqh),test=c_test,trainingLabel1 = trainingLabel,geneList = krasqh_300_genelist,out_dir_base = paste(basedir,paste("krasqh",300,"gene_list", sep="_"),sep='/'),single = 1)







testSig(sigProtein="kraswt",getGeneList=as.list(kraswt_100_genelist),trainingData = cbind(c_kras_gfp,c_kraswt), testData = c_test,geneList=NULL,trainingLabels = c(9,9))

foreach(numGenes = c(100,150,200,250,300)) %dopar% testSig(sigProtein="kraswt",numGenes = numGenes,trainingData = cbind(c_kras_gfp,c_kraswt), testData = c_test,geneList=NULL,trainingLabels = c(9,9))
foreach(numGenes = c(100,150,200,250,300)) %dopar% testSig(sigProtein="krasqh",numGenes = numGenes,trainingData = cbind(c_kras_gfp,c_krasqh), testData = c_test,geneList=NULL,trainingLabels = c(9,9))
foreach(numGenes = c(100,150,200,250,300)) %dopar% testSig(sigProtein="krasgv",numGenes = numGenes,trainingData = cbind(c_kras_gfp,c_krasgv), testData = c_test,geneList=NULL,trainingLabels = c(9,9))


luad_kras<-gatherFile("~/Dropbox/bild_signatures/kras_tcga_luad")
clinicals<-data.frame(fread("~/Desktop/PANCAN24/06_01_15_TCGA_24_548_Clinical_Variables_9264_Samples.txt"),check.names = F,row.names=T)
rownames(clinicals)
kras_mut_pancan_24<-t(clinicals["kras_mutation_found",3:ncol(clinicals)])
luad_preds_kras_mut<-merge_drop(kras_mut_pancan_24,luad_kras)
kras_yes<-subset(luad_preds_kras_mut,luad_preds_kras_mut$kras_mutation_found=="YES")
kras_no<-subset(luad_preds_kras_mut,luad_preds_kras_mut$kras_mutation_found=="NO")
pdf("~/Desktop/KRAS_validation.pdf")
for(i in 2:ncol(kras_yes)){
  #print(t.test(kras_yes[,i],kras_no[,i]))
  boxplot(kras_yes[,i],kras_no[,i],main=paste("p-value:",round(t.test(kras_yes[,i],kras_no[,i])$p.value,digits = 2),sep=" "),names=c("KRAS Mut", "KRAS WT"),xlab=colnames(kras_yes)[i])
}
dev.off()
# sam_kras_mut<-t(read.table("~/Dropbox/Datasets/SamKrasValidation/pancan_cleaned_filtered_annotated_KRAS_LungOnly.txt",header=1,row.names=1,sep='\t',check.names = F))
# dim(sam_kras_mut)
# sam_luad_kras<-data.frame(fread("~/Dropbox/Datasets/SamKrasValidation/LungAdenoRnaSeq.txt"),row.names=1,check.names = F)
# rownames(sam_luad_kras)

kras_mut<-t(read.table("~/Dropbox/bild_signatures/Datasets/kras_pancan12_mutations.txt",header=1,row.names=1,sep='\t',check.names = F))
sum(rownames(kras_mut)%in%rownames(sam_kras_mut))
shortnames=rownames(sam_kras_mut)
longnames=rownames(luad_kras)
new_longnames=NULL
counter=0
for (j in 1:length(longnames)){
  new_longnames[j]<-paste(strsplit(longnames[j],"-")[[1]][1:3],collapse = "-")
}
sum(new_longnames%in%rownames(kras_mut))
luad_kras$new_rownames=new_longnames
mut_pred_luad_new<-merge(luad_kras,kras_mut,by.x=19,by.y=0)
colnames(mut_pred_luad_new)=gsub(colnames(mut_pred_luad_new),pattern = "/pathway_activity_testset.csv/V1",replacement = "")
pdf("~/Dropbox/bild_signatures/kras_tcga_luad/kras_validation_boxplots.pdf")
for(i in 2:(ncol(mut_pred_luad_new)-1)){
  print(colnames(mut_pred_luad_new)[i])
  tmp=t.test(mut_pred_luad_new[,i]~mut_pred_luad_new[,ncol(mut_pred_luad_new)])
  print(tmp)
  boxplot(mut_pred_luad_new[,i]~mut_pred_luad_new[,ncol(mut_pred_luad_new)],main=paste(colnames(mut_pred_luad_new)[i],paste("p-val=",round(tmp$p.value,digits = 2),sep=" "),sep="\n"),names=c("KRAS_WT","KRAS_MUT"))
}
dev.off()




####
#source('~/Dropbox/bild_signatures/bild_signatures/Rmarkdowns_scripts//Key_ASSIGN_functions.Rmd', echo=TRUE)
setwd("~/Dropbox/bild_signatures/kras/")
filenames_icbp_multi<-system("ls */*/pathway_activity_testset*", intern=TRUE)
filenames_icbp_multi

for(i in 1:length(filenames_icbp_multi))
{
  f<-read.csv(filenames_icbp_multi[i], header=1,row.names=1) ###reading in the filess one at a time
  colnames(f)<-paste(filenames_icbp_multi[i],colnames(f),sep='/')
  if(i==1){
    data_icbp<-f
  }
  else{
    data_icbp<-cbind(data_icbp,f)
  }
}

head(data_icbp)
dim(data_icbp)
colnames(data_icbp)<-gsub(pattern = "/pathway_activity_testset.csv",replacement = "",x = colnames(data_icbp))
head(data_icbp)
rownames(data_icbp)[1:7]<-c("184A1","184B5","21MT1","21MT2","21NT","21PT","600MPE")
setwd("~/Dropbox/bild_signatures//Datasets")
drugs<-read.delim("ICBP_drugs.txt", header=1, sep='\t',row.names=1)
head(drugs);dim(drugs)
pred_drug<-merge_drop(data_icbp,drugs,by=0)
dim(pred_drug)
prot<-read.table("~/Dropbox/bild_signatures/bild_signatures/Datasets/proteomics.txt",sep='\t',header=1,row.names=1)
pred_prot<-merge_drop(data_icbp,prot)
cor_mat=p_mat=matrix(0,length(filenames_icbp_multi),70)
rownames(cor_mat)=rownames(p_mat)=colnames(pred_prot)[1:length(filenames_icbp_multi)]
colnames(cor_mat)=colnames(p_mat)=colnames(pred_prot)[(length(filenames_icbp_multi)+1):ncol(pred_prot)]

for(i in 1:length(filenames_icbp_multi)){
  for(j in 1:70){
    temp=cor.test(pred_prot[,i],pred_prot[,j+length(filenames_icbp_multi)],use="pairwise",method="spearman")
    print(j)
    print(temp)
    cor_mat[i,j]=temp$estimate
    p_mat[i,j]=temp$p.value
  }
}
colnames(p_mat)=paste(colnames(p_mat),"p_value",sep="_")
cor_p_mat<-cbind(cor_mat,p_mat)
order(colnames(cor_p_mat))
cor_p_mat<-cor_p_mat[,order(colnames(cor_p_mat))]
write.table(cor_p_mat,"~/Dropbox/bild_signatures/kras/ICBP_single_protein_kras_cor_p_mat_4_21.txt",col.names = NA,quote=F,sep='\t')
colnames(data_icbp)
data_icbp[,5]
##################
drugs<-read.delim("ICBP_drugs.txt", header=1, sep='\t',row.names=1)
icbp_drug<-merge_drop(data_icbp,drugs)
colnames(icbp_drug)
cor_mat=p_mat=matrix(0,length(filenames_icbp_multi),90)
rownames(cor_mat)=rownames(p_mat)=colnames(icbp_drug)[1:length(filenames_icbp_multi)]
colnames(cor_mat)=colnames(p_mat)=colnames(icbp_drug)[(length(filenames_icbp_multi)+11):ncol(icbp_drug)]

for(i in 1:length(filenames_icbp_multi)){
  for(j in 1:90){
    temp=cor.test(icbp_drug[,i],icbp_drug[,(j+length(filenames_icbp_multi)+10)],use="pairwise",method="spearman")
    print(j)
    print(temp)
    cor_mat[i,j]=temp$estimate
    p_mat[i,j]=temp$p.value
  }
}
#write.table(cor_mat,"~/Dropbox/bild_signatures/kras/ICBP_single_cor_drug_mat_4_21.txt",col.names = NA,quote=F,sep='\t')
colnames(p_mat)=paste(colnames(p_mat),"p_value",sep="_")
write.table(p_mat,"~/Dropbox/bild_signatures/kras/ICBP_single_p_drug_mat_4_21.txt",col.names = NA,quote=F,sep='\t')
cor_p_mat<-cbind(cor_mat,p_mat)
order(colnames(cor_p_mat))
cor_p_mat<-cor_p_mat[,order(colnames(cor_p_mat))]
write.table(cor_p_mat,"~/Dropbox/bild_signatures/kras/ICBP_single_cor_p_mat_6_18.txt",col.names = NA,quote=F,sep='\t')
#######################KRAS on TCGA LUAD#########
library(data.table)
gfp_kras<-read.table("~/Dropbox/Datasets/36_GFP_KRAS_TPMlog2.txt", sep='\t', header=1, row.names=1)
head(gfp_kras)
tcga<-data.frame(fread("~/Dropbox/bild_signatures/Datasets/PANCAN24_BRCA_1119_TPMlog2.txt"),check.names=F,row.names = 1)
sub=c(9,9,9,9,1119)
expr_f <-gfp_kras[apply(gfp_kras[,1:36]==0,1,mean) < 0.85,]
expr<-merge_drop(expr_f,tcga)
pcaplot(expr,sub)
bat<-as.matrix(cbind(colnames(expr),c(rep(1,ncol(gfp_kras)),rep(2,ncol(tcga)))))
combat_expr<-ComBat(dat=expr, batch=bat[,2], mod=NULL, numCovs=NULL)
pcaplot(combat_expr,sub)

# c_gfp<-subset(combat_expr, select=GFP.1:GFP.12)
# c_akt<-subset(combat_expr, select=AKT.1:AKT.6)
# c_bad<-subset(combat_expr, select=BAD.1:BAD.6)
# c_her2<-subset(combat_expr, select=HER2.1:HER2.6)
# c_igf1r<-subset(combat_expr, select=IGF1R.1:IGF1R.6)
# c_raf<-subset(combat_expr, select=RAF.1:RAF.6)
# c_erk<-subset(combat_expr, select=ERK.1:ERK.6)
# train_egfr<-combat_expr[,1:12]
c_kras_gfp<-subset(combat_expr,select=GFP.31:GFP.39)
c_kraswt<-subset(combat_expr,select=KRASWT.1:KRASWT.9)
c_krasqh<-subset(combat_expr,select=KRASQH.1:KRASQH.9)
c_krasgv<-subset(combat_expr,select=KRASGV.1:KRASGV.9)
c_test<-combat_expr[,(ncol(gfp_kras)+1):ncol(combat_expr)]
colnames(c_test)

#######getting the genelist#######
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/akt_75_gene_list/adapB_single/output.rda")
akt_75_genelist<-output.data$processed.data$diffGeneList  
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/single/akt_150_gene_list/adapB_single/output.rda")
akt_150_genelist<-output.data$processed.data$diffGeneList 
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/bad_200_gene_list/adapB_single/output.rda")
bad_200_genelist<-output.data$processed.data$diffGeneList  
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/igf1r_75_gene_list/adapB_single/output.rda")
igf1r_75_genelist<-output.data$processed.data$diffGeneList  
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/erk_250_gene_list/adapB_single/output.rda")
erk_250_genelist<-output.data$processed.data$diffGeneList  
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/single/erk_400_gene_list/adapB_single/output.rda")
erk_400_genelist<-output.data$processed.data$diffGeneList  
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/her2_15_gene_list/adapB_single/output.rda")
her2_15_genelist<-output.data$processed.data$diffGeneList  
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/egfr_25_gene_list/adapB_single/output.rda")
egfr_25_genelist<-output.data$processed.data$diffGeneList  
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/raf_100_gene_list/adapB_single/output.rda")
raf_100_genelist<-output.data$processed.data$diffGeneList


load("~/Dropbox/bild_signatures/kras/krasqh_300_gene_list/adapB_single/output.rda")
krasqh_300_genelist<-output.data$processed.data$diffGeneList
load("~/Dropbox/bild_signatures/kras/krasgv_300_gene_list/adapB_single/output.rda")
krasgv_300_genelist<-output.data$processed.data$diffGeneList
load("~/Dropbox/bild_signatures/kras/kraswt_300_gene_list/adapB_single/output.rda")
kraswt_300_genelist<-output.data$processed.data$diffGeneList




#single_pathway_best<-data_icbp[,c("akt_75_gene_list/adap_adap_single/akt","her2_15_gene_list/adap_adap_single/her2","igf1r_75_gene_list/adap_adap_single/igf1r","bad_200_gene_list/adap_adap_single/bad","erk_250_gene_list/adap_adap_single/erk","raf_100_gene_list/adap_adap_single/raf","egfr_25_gene_list/adap_adap_single/egfr")]
#best_pathway_predictions<-pred_drug[,c("akt_75_gene_list/adap_adap_single/akt","her2_25_gene_list/adap_adap_single/her2","igf1r_75_gene_list/adap_adap_single/igf1r","bad_200_gene_list/adap_adap_single/bad","erk_250_gene_list/adap_adap_single/erk","raf_100_gene_list/adap_adap_single/raf","egfr_25_gene_list/adap_adap_single/egfr")]
#single_pathway_best<-c("akt_75_gene_list/adap_adap_single/akt","bad_200_gene_list/adap_adap_single/bad","her2_15_gene_list/adap_adap_single/her2","igf1r_75_gene_list/adap_adap_single/igf1r","erk_200_gene_list/adap_adap_single/erk","raf_100_gene_list/adap_adap_single/raf","egfr_25_gene_list/adap_adap_single/egfr")
#multi_pathway_best<-c("akt_bad_her2_raf/adap_adap_multi/akt","akt_bad_her2_erk/adap_adap_multi/bad","her2_igf1r/adap_adap_multi/her2","igf1r_raf_her2_erk/adap_adap_multi/igf1r","erk_bad/adap_adap_multi/erk","akt_bad_igf1r_raf/adapB_multi/raf")

#####validation in TCGA LUAD 541 samples based on KRAS mutation status. Hypothesis is, mutant sample will have higher kras activity than wild type kras samples;
###KRASWT signature was used as control and expected to see no difference.
basedir="~/Dropbox/bild_signatures/kras_tcga_brca/bioconductor"
dir.create(basedir)
testSig <- function(sigProtein, numGenes=NA, geneList =NULL, trainingData, testData, trainingLabels)
{
  names(sigProtein)<-sigProtein
  trainingLabel<-list(control=list(sigProtein=1:trainingLabels[1]),sigProtein=(trainingLabels[1]+1):(trainingLabels[1]+trainingLabels[2]))
  names(trainingLabel$control)=names(trainingLabel)[2]=names(sigProtein)
  if(is.na(numGenes)){
    sub_dir<-paste(basedir,paste(sigProtein,"gene_list", sep="_"),sep='/')
  }
  else{
    sub_dir<-paste(basedir,paste(sigProtein,numGenes,"gene_list", sep="_"),sep='/')
  }
  dir.create(sub_dir)
  assign_easy_multi(trainingData = trainingData,test=testData,trainingLabel1 = trainingLabel,g=numGenes,geneList = geneList,out_dir_base = sub_dir,single = 1)
}

foreach(numGenes = c(300)) %dopar% testSig(sigProtein="kraswt",numGenes = numGenes,trainingData = cbind(c_kras_gfp,c_kraswt), testData = c_test,geneList=NULL,trainingLabels = c(9,9))
foreach(numGenes = c(300)) %dopar% testSig(sigProtein="krasqh",numGenes = numGenes,trainingData = cbind(c_kras_gfp,c_krasqh), testData = c_test,geneList=NULL,trainingLabels = c(9,9))
foreach(numGenes = c(300)) %dopar% testSig(sigProtein="krasgv",numGenes = numGenes,trainingData = cbind(c_kras_gfp,c_krasgv), testData = c_test,geneList=NULL,trainingLabels = c(9,9))
