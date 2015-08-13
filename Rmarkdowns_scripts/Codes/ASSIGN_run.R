args=(commandArgs(TRUE))
testFile = commandArgs()[7]
basedir= commandArgs()[8]
dir.create(basedir)
library(devtools)
install_github("wevanjohnson/ASSIGN",ref="adapt_gene_only_ver2")
library(ASSIGN)
library(foreach)
library(sva)
library(data.table)
library(doParallel)
library(doMC)  
registerDoMC(8)
#library(BatchQC)
#cl<-makeCluster(2)
#registerDoParallel(cl)


source("Key_ASSIGN_functions.Rmd")

setwd("/data2/Moom_datasets/")
expr<-as.matrix(read.table("GFP18_AKT_BAD_HER2_IGF1R_RAF_ERK.tpmlog",sep='\t',row.names=1,header=1))
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
control_egfr_l<-read.table("18_GFP_EGFR_TPMlog2.txt", sep='\t', header=1, row.names=1)
gfp_egfr_multi_f <- merge_drop(control_egfr_l,expr_all_f)
dim(gfp_egfr_multi_f)

#gfp_kras<-read.table("36_GFP_KRAS_TPMlog2.txt", sep='\t', header=1, row.names=1)
#head(gfp_kras)
#gfp_egfr_kras_multi_f<-merge_drop(gfp_egfr_multi_f,gfp_kras)
#dim(gfp_egfr_kras_multi_f)
test<-data.frame(fread(testFile), check.names=F,row.names=1)
dim(test)
expr_all_test_f<-merge_drop(gfp_egfr_multi_f,test,by=0)

dim(expr_all_test_f)
colnames(expr_all_test_f)
sub<-c(6,6,12,6,6,5,6,6,6,ncol(test))

pdf(file=paste(basedir,'pcaPlots.pdf',sep='_'))
pcaplot(expr_all_test_f,sub)
bat1<-as.matrix(cbind(colnames(expr_all_test_f),c(rep(1,ncol(control_egfr_l)),rep(2,ncol(expr_all_f)),rep(3,ncol(test)))))
#bat1
combat_expr1<-ComBat(dat=expr_all_test_f, batch=bat1[,2], mod=NULL, numCovs=NULL)
pcaplot(combat_expr1,sub)
dev.off()
c_gfp<-subset(combat_expr1, select=GFP.1:GFP.12)
c_akt<-subset(combat_expr1, select=AKT.1:AKT.6)
c_bad<-subset(combat_expr1, select=BAD.1:BAD.6)
c_her2<-subset(combat_expr1, select=HER2.1:HER2.6)
c_igf1r<-subset(combat_expr1, select=IGF1R.1:IGF1R.6)
c_raf<-subset(combat_expr1, select=RAF.1:RAF.6)
c_erk<-subset(combat_expr1, select=ERK.1:ERK.6)
train_egfr<-combat_expr1[,1:12]
c_test<-combat_expr1[,(ncol(gfp_egfr_multi_f)+1):ncol(combat_expr1)]
colnames(c_test)

#######getting the genelist#######
load("/data2/Moom_datasets/GeneLists/akt_75_gene_list/adapB_single/output.rda")
akt_75_genelist<-output.data$processed.data$diffGeneList  
load("/data2/Moom_datasets/GeneLists/bad_200_gene_list/adapB_single/output.rda")
bad_200_genelist<-output.data$processed.data$diffGeneList  
load("/data2/Moom_datasets/GeneLists/igf1r_75_gene_list/adapB_single/output.rda")
igf1r_75_genelist<-output.data$processed.data$diffGeneList  
load("/data2/Moom_datasets/GeneLists/erk_250_gene_list/adapB_single/output.rda")
erk_250_genelist<-output.data$processed.data$diffGeneList  
load("/data2/Moom_datasets/GeneLists/her2_15_gene_list/adapB_single/output.rda")
her2_15_genelist<-output.data$processed.data$diffGeneList  
load("/data2/Moom_datasets/GeneLists/egfr_25_gene_list/adapB_single/output.rda")
egfr_25_genelist<-output.data$processed.data$diffGeneList  
load("/data2/Moom_datasets/GeneLists/raf_100_gene_list/adapB_single/output.rda")
raf_100_genelist<-output.data$processed.data$diffGeneList


#foreach(sigProtein = c('akt_75')) %dopar% testSig(sigProtein=sigProtein,numGenes=NA,trainingData = cbind(c_gfp,eval(parse(text=paste('c',strsplit(sigProtein,"_")[[1]][1],sep="_")))), testData = c_test,geneList=eval(parse(text=paste(sigProtein,"genelist",sep="_"))),trainingLabels = c(12,6))
foreach(sigProtein = c('akt_75','bad_200','raf_100','igf1r_75','erk_250')) %dopar% testSig(sigProtein=sigProtein,numGenes=NA,trainingData = cbind(c_gfp,eval(parse(text=paste('c',strsplit(sigProtein,"_")[[1]][1],sep="_")))), testData = c_test,geneList=eval(parse(text=paste(sigProtein,"genelist",sep="_"))),trainingLabels = c(12,6))
testSig("her2_15", geneList=eval(parse(text=paste("her2_15","genelist",sep="_"))), trainingData = cbind(c_gfp,eval(parse(text=paste('c',"her2",sep="_")))), testData = c_test, trainingLabels = c(12,5))
testSig("egfr_25", geneList=eval(parse(text=paste("egfr_25","genelist",sep="_"))), trainingData = train_egfr, testData = c_test, trainingLabels = c(6,6))

##best_her2
trainingLabel<-list(control=list(her2=1:12,igf1r=1:12),her2=13:17, igf1r=18:23)
her2_igf1r_genelist<-c(eval(parse(text=paste("her2_15","genelist",sep="_"))),eval(parse(text=paste("igf1r_75","genelist",sep="_"))))
her2_igf1r_genelist<-c(her2_15_genelist,igf1r_75_genelist)
names(her2_igf1r_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="her2_igf1r", geneList=her2_igf1r_genelist, trainingData = cbind(c_gfp,c_her2,c_igf1r), testData = c_test, trainingLabel = trainingLabel)

##best_akt
trainingLabel<-list(control=list(akt=1:12,bad=1:12),akt=13:18, bad=19:24)
akt_bad_genelist<-c(akt_75_genelist,bad_200_genelist)
names(akt_bad_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="akt_bad", geneList=akt_bad_genelist, trainingData = cbind(c_gfp,c_akt,c_bad), testData = c_test, trainingLabel = trainingLabel)


##best_bad: akt_bad_her2_erk
trainingLabel<-list(control=list(akt=1:12,bad=1:12,her2=1:12,erk=1:12),akt=13:18, bad=19:24,her2=25:29,erk=30:35)
akt_bad_her2_erk_genelist<-c(akt_75_genelist,bad_200_genelist,her2_15_genelist,erk_250_genelist)
names(akt_bad_her2_erk_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="akt_bad_her2_erk", geneList=akt_bad_her2_erk_genelist, trainingData = cbind(c_gfp,c_akt,c_bad,c_her2,c_erk), testData = c_test, trainingLabel = trainingLabel)

##best_raf:akt_bad_igf1r_raf
trainingLabel<-list(control=list(akt=1:12,bad=1:12,igf1r=1:12,raf=1:12),akt=13:18, bad=19:24,igf1r=25:30,raf=31:36)
akt_bad_igf1r_raf_genelist<-c(akt_75_genelist,bad_200_genelist,igf1r_75_genelist,raf_100_genelist)
names(akt_bad_igf1r_raf_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="akt_bad_igf1r_raf", geneList=akt_bad_igf1r_raf_genelist, trainingData = cbind(c_gfp,c_akt,c_bad,c_igf1r,c_raf), testData = c_test, trainingLabel = trainingLabel)

##best_igf1r:igf1r_raf_her2_erk
trainingLabel<-list(control=list(igf1r=1:12,raf=1:12,her2=1:12,erk=1:12),igf1r=13:18, raf=19:24,her2=25:29,erk=30:35)
igf1r_raf_her2_erk_genelist<-c(igf1r_75_genelist,raf_100_genelist,her2_15_genelist,erk_250_genelist)
names(igf1r_raf_her2_erk_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="igf1r_raf_her2_erk", geneList=igf1r_raf_her2_erk_genelist, trainingData = cbind(c_gfp,c_igf1r,c_raf,c_her2,c_erk), testData = c_test, trainingLabel = trainingLabel)

##best_erk:erk_bad
trainingLabel<-list(control=list(erk=1:12,bad=1:12),erk=13:18, bad=19:24)
erk_bad_genelist<-c(erk_250_genelist,bad_200_genelist)
names(erk_bad_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="erk_bad", geneList=erk_bad_genelist, trainingData = cbind(c_gfp,c_erk,c_bad), testData = c_test, trainingLabel = trainingLabel)

