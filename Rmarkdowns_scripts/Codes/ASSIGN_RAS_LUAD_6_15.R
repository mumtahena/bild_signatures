library(devtools)
install_github("wevanjohnson/ASSIGN")
library(ASSIGN)
library(foreach)
library(sva)
library(data.table)
library(doParallel)

cl<-makeCluster(2)
registerDoParallel(cl)
source("Key_ASSIGN_functions.Rmd")
setwd("/data2/Moom_datasets/")
gfp_kras<-read.table("36_GFP_KRAS_TPMlog2.txt", sep='\t', header=1, row.names=1)
head(gfp_kras)
expr_all_f <-gfp_kras[apply(gfp_kras==0,1,mean) < 0.80,]
test<-data.frame(fread("/data2/Moom_datasets/PANCAN24/PANCAN24_LUAD_541_TPMlog2.txt"), check.names=F,row.names=1)
dim(test)
colnames(test)
expr_all_test_f<-merge_drop(expr_all_f,test,by=0)
sub<-c(9,9,9,9,ncol(test))
pdf("~/kras_LUAD.pdf")
pcaplot(expr_all_test_f,sub)#,scale=F,center=F)
bat1<-as.matrix(cbind(colnames(expr_all_test_f),c(rep(1,ncol(gfp_kras)),rep(2,ncol(test)))))
bat1
combat_expr<-ComBat(dat=expr_all_test_f, batch=bat1[,2], mod=NULL, numCovs=NULL)
pcaplot(combat_expr,sub)#,scale=F, center=F)
dev.off()


c_kras_gfp<-subset(combat_expr,select=GFP.31:GFP.39)
c_kraswt<-subset(combat_expr,select=KRASWT.1:KRASWT.9)
c_krasqh<-subset(combat_expr,select=KRASQH.1:KRASQH.9)
c_krasgv<-subset(combat_expr,select=KRASGV.1:KRASGV.9)
c_test<-combat_expr[,(ncol(gfp_kras)+1):ncol(combat_expr)]
basedir="~/Kras_pred_LUAD"
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
#library(doParallel)
#registerDoParallel(cores=2)
foreach(numGenes = c(100,150)) %dopar% testSig(sigProtein="kraswt",numGenes = numGenes,trainingData = cbind(c_kras_gfp,

#foreach(numGenes = c(100,150,200,250,300)) %dopar% testSig(sigProtein="kraswt",numGenes = numGenes,trainingData = cbind(c_kras_gfp,c_kraswt), testData = c_test,geneList=NULL,trainingLabels = c(9,9))
#foreach(numGenes = c(100,150,200,250,300)) %dopar% testSig(sigProtein="krasqh",numGenes = numGenes,trainingData = cbind(c_kras_gfp,c_krasqh), testData = c_test,geneList=NULL,trainingLabels = c(9,9))
#foreach(numGenes = c(100,150,200,250,300)) %dopar% testSig(sigProtein="krasgv",numGenes = numGenes,trainingData = cbind(c_kras_gfp,c_krasgv), testData = c_test,geneList=NULL,trainingLabels = c(9,9))




stopCluster(cl)
