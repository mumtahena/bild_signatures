library(sva)
library(data.table)
library(devtools)
install_github("wevanjohnson/ASSIGN",ref="adapt_gene_only_ver2")
library(ASSIGN)
#install_github("mani2012/BatchQC", build_vignettes=TRUE)
#library(BatchQC)

source("~/Dropbox/bild_signatures//bild_signatures/Rmarkdowns_scripts/Key_ASSIGN_functions.Rmd")
#setwd("~/Documents/ThesisWork/GitRepos/bild_signature_validation_old_repo/Datasets")
setwd("/Users/mumtahenarahman/Dropbox/bild_signature/Datasets/")
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
multi<-read.table("~/Dropbox/egfr_opt/hmec_multi_egfr_mek_FeatureCount.tpmlog", sep='\t', header=1, row.names=1)
head(multi)
dim(multi)
multi<-multi[,order(colnames(multi))]
head(multi)
colnames(multi)
control_l<-multi[,1:6]
egfr_l<-multi[,7:12]
control_egfr_l<-cbind(control_l,egfr_l)
gfp_egfr_multi_f <- merge_drop(control_egfr_l,expr_all_f)
dim(gfp_egfr_multi_f)

icbp<-as.matrix(read.table("~/Dropbox/bild_signatures/Datasets/icbp_Rsubread_tpmlog.txt", sep='\t', stringsAsFactors=FALSE, header=1, row.names=1))
gfp_egfr_icbp_f<-merge_drop(gfp_egfr_multi_f,icbp,by=0)
colnames(gfp_egfr_icbp_f)
sub<-c(6,6,12,6,6,5,6,6,6,55)
pdf("pca_plots_all_7_pathways_5_9.pdf")
pcaplot(mat = gfp_egfr_icbp_f,sub = sub)
length(colnames(gfp_egfr_icbp_f))
bat1<-as.matrix(cbind(c(colnames(gfp_egfr_multi_f),colnames(icbp)),as.numeric(c(rep(1,12),rep(2,47),rep(3,55)))))
conditions<-c(rep("Control",6),rep("EGFR",6),rep("Control",12),rep("AKT",6),rep("BAD",6),rep("HER2",5),rep("IGF1R",6),rep("RAF",6),rep("ERK",7),rep("ICBP",55))
modmatrix = model.matrix(~as.factor(conditions), data=gfp_egfr_multi_f)
combat_expr1<-ComBat(dat=gfp_egfr_icbp_f, batch=bat1[,2], mod=NULL, numCovs=NULL)
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
c_test<-combat_expr1[,60:114]
colnames(c_test)

#######getting the genelist#######
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

########running multipathway genelist based predictions with ASSIGN#######
pathways<-c("akt","bad","her2","igf1r","erk","egfr","raf")
pathways_2<-combn(pathways,2)
egfr_combos_2=NULL

for( i in 1:ncol(pathways_2)){
  if("egfr"%in%pathways_2[,i]){
    if(is.null(egfr_combos_2)){
      egfr_combos_2=pathways_2[,i]
    }
    else{
      egfr_combos_2=cbind(egfr_combos_2,pathways_2[,i])
    }
  }
}

foreach(sigProtein = c('akt_75','akt_150','bad_200','raf_100','igf1r_75','erk_250','erk_400')) %dopar% testSig(sigProtein,trainingData = cbind(c_gfp,eval(parse(text=paste('c',strsplit(sigProtein,"_")[[1]][1],sep="_")))), testData = c_test,geneList=eval(parse(text=paste(sigProtein,"genelist",sep="_"))),trainingLabels = c(12,6))
testSig("her2_15", geneList=eval(parse(text=paste("her2_15","genelist",sep="_"))), trainingData = cbind(c_gfp,eval(parse(text=paste('c',"her2",sep="_")))), testData = c_test, trainingLabels = c(12,5))
testSig("egfr_25", geneList=eval(parse(text=paste("egfr_25","genelist",sep="_"))), trainingData = train_egfr, testData = c_test, trainingLabels = c(6,6))
##multipathway
multi_pathway_best<-c("akt_bad/adap_adap_multi/akt","akt_bad_her2_erk/adap_adap_multi/bad","her2_igf1r/adap_adap_multi/her2","igf1r_raf_her2_erk/adap_adap_multi/igf1r","erk_bad/adap_adap_multi/erk","akt_bad_igf1r_raf/adapB_multi/raf")


testSig_multi <- function(sigProteins, numGenes=NA,geneList =NULL, trainingData, testData, trainingLabel)
{
  sub_dir<-paste(basedir,paste(sigProteins,"gene_list", sep="_"),sep='/')
  dir.create(sub_dir)
  assign_easy_multi(trainingData = trainingData,test=testData,trainingLabel1 = trainingLabel,g=numGenes,geneList = geneList,out_dir_base = sub_dir)
}

##akt_egfr
basedir="~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/"
trainingLabel<-list(control=list(akt=1:12,egfr=19:24),akt=13:18, egfr=25:30)
akt_egfr_genelist<-c(eval(parse(text=paste("akt_75","genelist",sep="_"))),eval(parse(text=paste("egfr_25","genelist",sep="_"))))
names(akt_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="akt_egfr", geneList=akt_egfr_genelist, trainingData = cbind(c_gfp,c_akt,train_egfr), testData = c_test, trainingLabel = trainingLabel)

##bad_egfr
trainingLabel<-list(control=list(bad=1:12,egfr=19:24),bad=13:18, egfr=25:30)
bad_egfr_genelist<-c(eval(parse(text=paste("bad_200","genelist",sep="_"))),eval(parse(text=paste("egfr_25","genelist",sep="_"))))
names(bad_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="bad_egfr", geneList=bad_egfr_genelist, trainingData = cbind(c_gfp,c_bad,train_egfr), testData = c_test, trainingLabel = trainingLabel)

##her2_egfr
trainingLabel<-list(control=list(her2=1:12,egfr=18:23),her2=13:17, egfr=24:29)
her2_egfr_genelist<-c(eval(parse(text=paste("her2_15","genelist",sep="_"))),eval(parse(text=paste("egfr_25","genelist",sep="_"))))
names(her2_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="her2_egfr", geneList=her2_egfr_genelist, trainingData = cbind(c_gfp,c_her2,train_egfr), testData = c_test, trainingLabel = trainingLabel)

##igf1r_egfr
trainingLabel<-list(control=list(igf1r=1:12,egfr=19:24),akt=13:18, egfr=25:30)
igf1r_egfr_genelist<-c(eval(parse(text=paste("igf1r_75","genelist",sep="_"))),eval(parse(text=paste("egfr_25","genelist",sep="_"))))
names(igf1r_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="igf1r_egfr", geneList=igf1r_egfr_genelist, trainingData = cbind(c_gfp,c_igf1r,train_egfr), testData = c_test, trainingLabel = trainingLabel)

##erk_egfr
trainingLabel<-list(control=list(erk=1:12,egfr=19:24),erk=13:18, egfr=25:30)
erk_egfr_genelist<-c(eval(parse(text=paste("erk_250","genelist",sep="_"))),eval(parse(text=paste("egfr_25","genelist",sep="_"))))
names(erk_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="erk_egfr", geneList=erk_egfr_genelist, trainingData = cbind(c_gfp,c_erk,train_egfr), testData = c_test, trainingLabel = trainingLabel)

## raf_egfr
trainingLabel<-list(control=list(raf=1:12,egfr=19:24),raf=13:18, egfr=25:30)
raf_egfr_genelist<-c(eval(parse(text=paste("raf_100","genelist",sep="_"))),eval(parse(text=paste("egfr_25","genelist",sep="_"))))
names(raf_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="raf_egfr", geneList=bad_egfr_genelist, trainingData = cbind(c_gfp,c_raf,train_egfr), testData = c_test, trainingLabel = trainingLabel)


##### [,1]  [,2]   [,3]    [,4]  [,5]   [,6]  [,7]   [,8]    [,9]  [,10]  [,11] [,12]  
#[1,] "AKT" "AKT"  "AKT"   "AKT" "AKT"  "AKT" "BAD"  "BAD"   "BAD" "BAD"  "BAD" "HER2" 
#[2,] "BAD" "HER2" "IGF1R" "ERK" "EGFR" "RAF" "HER2" "IGF1R" "ERK" "EGFR" "RAF" "IGF1R"
#####[,13]  [,14]  [,15]  [,16]   [,17]   [,18]   [,19]  [,20] [,21] 
#[1,] "HER2" "HER2" "HER2" "IGF1R" "IGF1R" "IGF1R" "ERK"  "ERK" "EGFR"
#[2,] "ERK"  "EGFR" "RAF"  "ERK"   "EGFR"  "RAF"   "EGFR" "RAF" "RAF" 
pathways_3<-combn(pathways,3)

#1. HER2 & AKT & EGFR

trainingLabel<-list(control=list(akt=1:12,her2=1:12,egfr=24:29),akt=13:18,her2=19:23,egfr=30:35)
akt_her2_egfr_genelist<-c(akt_75_genelist,her2_15_genelist,egfr_25_genelist)
names(akt_her2_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="akt_her2_egfr", geneList=akt_her2_egfr_genelist, trainingData = cbind(c_gfp,c_akt,c_her2,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#2. HER2 & BAD & EGFR
trainingLabel<-list(control=list(bad=1:12,her2=1:12,egfr=24:29),bad=13:18,her2=19:23,egfr=30:35)
bad_her2_egfr_genelist<-c(bad_200_genelist,her2_15_genelist,egfr_25_genelist)
names(bad_her2_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="bad_her2_egfr", geneList=bad_her2_egfr_genelist, trainingData = cbind(c_gfp,c_bad,c_her2,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#3. HER2 & IGF1R & EGFR
trainingLabel<-list(control=list(igf1r=1:12,her2=1:12,egfr=24:29),igf1r=13:18,her2=19:23,egfr=30:35)
igf1r_her2_egfr_genelist<-c(igf1r_75_genelist,her2_15_genelist,egfr_25_genelist)
names(igf1r_her2_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="igf1r_her2_egfr", geneList=bad_her2_egfr_genelist, trainingData = cbind(c_gfp,c_igf1r,c_her2,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#4. HER2 & ERK & EGFR
trainingLabel<-list(control=list(erk=1:12,her2=1:12,egfr=24:29),erk=13:18,her2=19:23,egfr=30:35)
erk_her2_egfr_genelist<-c(igf1r_75_genelist,her2_15_genelist,egfr_25_genelist)
names(erk_her2_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="erk_her2_egfr", geneList=erk_her2_egfr_genelist, trainingData = cbind(c_gfp,c_erk,c_her2,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#5. HER2 & RAF & EGFR
trainingLabel<-list(control=list(raf=1:12,her2=1:12,egfr=24:29),raf=13:18,her2=19:23,egfr=30:35)
raf_her2_egfr_genelist<-c(raf_100_genelist,her2_15_genelist,egfr_25_genelist)
names(raf_her2_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="raf_her2_egfr", geneList=raf_her2_egfr_genelist, trainingData = cbind(c_gfp,c_raf,c_her2,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#6. AKT & BAD & EGFR
trainingLabel<-list(control=list(akt=1:12,bad=1:12,egfr=25:30),akt=13:18,bad=19:24,egfr=31:36)
akt_bad_egfr_genelist<-c(akt_75_genelist,bad_200_genelist,egfr_25_genelist)
names(akt_bad_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="akt_bad_egfr", geneList=akt_bad_egfr_genelist, trainingData = cbind(c_gfp,c_akt,c_bad,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#7. AKT & IGF1R & EGFR
trainingLabel<-list(control=list(akt=1:12,igf1r=1:12,egfr=25:30),akt=13:18,igf1r=19:24,egfr=31:36)
akt_igf1r_egfr_genelist<-c(akt_75_genelist,igf1r_75_genelist,egfr_25_genelist)
names(akt_igf1r_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="akt_igf1r_egfr", geneList=akt_igf1r_egfr_genelist, trainingData = cbind(c_gfp,c_akt,c_igf1r,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#8. AKT & ERK & EGFR
trainingLabel<-list(control=list(akt=1:12,erk=1:12,egfr=25:30),akt=13:18,erk=19:24,egfr=31:36)
akt_erk_egfr_genelist<-c(akt_75_genelist,erk_250_genelist,egfr_25_genelist)
names(akt_erk_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="akt_erk_egfr", geneList=akt_erk_egfr_genelist, trainingData = cbind(c_gfp,c_akt,c_erk,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#9. AKT & RAF & EGFR
trainingLabel<-list(control=list(akt=1:12,raf=1:12,egfr=25:30),akt=13:18,raf=19:24,egfr=31:36)
akt_raf_egfr_genelist<-c(akt_75_genelist,raf_100_genelist,egfr_25_genelist)
names(akt_raf_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="akt_raf_egfr", geneList=akt_raf_egfr_genelist, trainingData = cbind(c_gfp,c_akt,c_raf,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#10. BAD & IGF1R & EGFR
trainingLabel<-list(control=list(bad=1:12,igf1r=1:12,egfr=25:30),bad=13:18,igf1r=19:24,egfr=31:36)
bad_igf1r_egfr_genelist<-c(bad_200_genelist,igf1r_75_genelist,egfr_25_genelist)
names(bad_igf1r_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="bad_igf1r_egfr", geneList=bad_igf1r_egfr_genelist, trainingData = cbind(c_gfp,c_bad,c_igf1r,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#11. BAD & ERK & EGFR
trainingLabel<-list(control=list(bad=1:12,erk=1:12,egfr=25:30),bad=13:18,erk=19:24,egfr=31:36)
bad_erk_egfr_genelist<-c(bad_200_genelist,erk_250_genelist,egfr_25_genelist)
names(bad_erk_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="bad_erk_egfr", geneList=bad_erk_egfr_genelist, trainingData = cbind(c_gfp,c_bad,c_erk,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#12. BAD & RAF & EGFR
trainingLabel<-list(control=list(bad=1:12,raf=1:12,egfr=25:30),bad=13:18,raf=19:24,egfr=31:36)
bad_raf_egfr_genelist<-c(bad_200_genelist,raf_100_genelist,egfr_25_genelist)
names(bad_raf_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="bad_raf_egfr", geneList=bad_raf_egfr_genelist, trainingData = cbind(c_gfp,c_bad,c_raf,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#13. IGF1R & RAF & EGFR
trainingLabel<-list(control=list(igf1r=1:12,raf=1:12,egfr=25:30),igf1r=13:18,raf=19:24,egfr=31:36)
igf1r_raf_egfr_genelist<-c(igf1r_75_genelist,raf_100_genelist,egfr_25_genelist)
names(igf1r_raf_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="igf1r_raf_egfr", geneList=igf1r_raf_egfr_genelist, trainingData = cbind(c_gfp,c_igf1r,c_raf,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#14. IGF1R & ERK & EGFR
trainingLabel<-list(control=list(igf1r=1:12,erk=1:12,egfr=25:30),igf1r=13:18,erk=19:24,egfr=31:36)
igf1r_erk_egfr_genelist<-c(igf1r_75_genelist,erk_250_genelist,egfr_25_genelist)
names(igf1r_erk_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="igf1r_erk_egfr", geneList=igf1r_erk_egfr_genelist, trainingData = cbind(c_gfp,c_igf1r,c_erk,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#15. RAF & ERK & EGFR
trainingLabel<-list(control=list(raf=1:12,erk=1:12,egfr=25:30),raf=13:18,erk=19:24,egfr=31:36)
raf_erk_egfr_genelist<-c(raf_100_genelist,erk_250_genelist,egfr_25_genelist)
names(raf_erk_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="raf_erk_egfr", geneList=raf_erk_egfr_genelist, trainingData = cbind(c_gfp,c_raf,c_erk,train_egfr), testData = c_test, trainingLabel = trainingLabel)






######[,1]   [,2]    [,3]  [,4]   [,5]  [,6]    [,7]   [,8]   [,9]   [,10]   [,11]  
#[1,] "AKT"  "AKT"   "AKT" "AKT"  "AKT" "AKT"   "AKT"  "AKT"  "AKT"  "AKT"   "AKT"  
#[2,] "BAD"  "BAD"   "BAD" "BAD"  "BAD" "HER2"  "HER2" "HER2" "HER2" "IGF1R" "IGF1R"
#[3,] "HER2" "IGF1R" "ERK" "EGFR" "RAF" "IGF1R" "ERK"  "EGFR" "RAF"  "ERK"   "EGFR" 
#####[,12]   [,13]  [,14] [,15]  [,16]   [,17]  [,18]  [,19]  [,20]   [,21]   [,22]  
#[1,] "AKT"   "AKT"  "AKT" "AKT"  "BAD"   "BAD"  "BAD"  "BAD"  "BAD"   "BAD"   "BAD"  
#[2,] "IGF1R" "ERK"  "ERK" "EGFR" "HER2"  "HER2" "HER2" "HER2" "IGF1R" "IGF1R" "IGF1R"
#[3,] "RAF"   "EGFR" "RAF" "RAF"  "IGF1R" "ERK"  "EGFR" "RAF"  "ERK"   "EGFR"  "RAF"  
#####[,23]  [,24] [,25]  [,26]   [,27]   [,28]   [,29]  [,30]  [,31]  [,32]   [,33]  
#[1,] "BAD"  "BAD" "BAD"  "HER2"  "HER2"  "HER2"  "HER2" "HER2" "HER2" "IGF1R" "IGF1R"
#[2,] "ERK"  "ERK" "EGFR" "IGF1R" "IGF1R" "IGF1R" "ERK"  "ERK"  "EGFR" "ERK"   "ERK"  
#[3,] "EGFR" "RAF" "RAF"  "ERK"   "EGFR"  "RAF"   "EGFR" "RAF"  "RAF"  "EGFR"  "RAF"  
#####[,34]   [,35] 
#[1,] "IGF1R" "ERK" 
#[2,] "EGFR"  "EGFR"
#[3,] "RAF"   "RAF" 
pathways_4<-combn(pathways,4)
######[,1]    [,2]   [,3]   [,4]   [,5]    [,6]    [,7]    [,8]   [,9]  [,10]  [,11]  
#[1,] "AKT"   "AKT"  "AKT"  "AKT"  "AKT"   "AKT"   "AKT"   "AKT"  "AKT" "AKT"  "AKT"  
#[2,] "BAD"   "BAD"  "BAD"  "BAD"  "BAD"   "BAD"   "BAD"   "BAD"  "BAD" "BAD"  "HER2" 
#[3,] "HER2"  "HER2" "HER2" "HER2" "IGF1R" "IGF1R" "IGF1R" "ERK"  "ERK" "EGFR" "IGF1R"
#[4,] "IGF1R" "ERK"  "EGFR" "RAF"  "ERK"   "EGFR"  "RAF"   "EGFR" "RAF" "RAF"  "ERK"  

#1. HER2 & AKT & BAD & EGFR

trainingLabel<-list(control=list(akt=1:12,bad=1:12,her2=1:12,egfr=30:35),akt=13:18,bad=19:24,her2=25:29,egfr=36:41)
akt_bad_her2_egfr_genelist<-c(akt_75_genelist,bad_200_genelist,her2_15_genelist,egfr_25_genelist)
names(akt_bad_her2_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="akt_bad_her2_egfr", geneList=akt_bad_her2_egfr_genelist, trainingData = cbind(c_gfp,c_akt,c_bad,c_her2,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#2. HER2 & AKT & IGF1R & EGFR
trainingLabel<-list(control=list(akt=1:12,igf1r=1:12,her2=1:12,egfr=30:35),akt=13:18,igf1r=19:24,her2=25:29,egfr=36:41)
akt_igf1r_her2_egfr_genelist<-c(akt_75_genelist,igf1r_75_genelist,her2_15_genelist,egfr_25_genelist)
names(akt_igf1r_her2_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="akt_igf1r_her2_egfr", geneList=akt_igf1r_her2_egfr_genelist, trainingData = cbind(c_gfp,c_akt,c_igf1r,c_her2,train_egfr), testData = c_test, trainingLabel = trainingLabel)


##3. HER2 & AKT & RAF & EGFR
trainingLabel<-list(control=list(akt=1:12,raf=1:12,her2=1:12,egfr=30:35),akt=13:18,raf=19:24,her2=25:29,egfr=36:41)
akt_raf_her2_egfr_genelist<-c(akt_75_genelist,raf_100_genelist,her2_15_genelist,egfr_25_genelist)
names(akt_raf_her2_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="akt_raf_her2_egfr", geneList=akt_raf_her2_egfr_genelist, trainingData = cbind(c_gfp,c_akt,c_raf,c_her2,train_egfr), testData = c_test, trainingLabel = trainingLabel)

##4. HER2 & AKT & ERK & EGFR
trainingLabel<-list(control=list(akt=1:12,erk=1:12,her2=1:12,egfr=30:35),akt=13:18,erk=19:24,her2=25:29,egfr=36:41)
akt_erk_her2_egfr_genelist<-c(akt_75_genelist,erk_250_genelist,her2_15_genelist,egfr_25_genelist)
names(akt_erk_her2_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="akt_erk_her2_egfr", geneList=akt_erk_her2_egfr_genelist, trainingData = cbind(c_gfp,c_akt,c_erk,c_her2,train_egfr), testData = c_test, trainingLabel = trainingLabel)


#5. HER2 & BAD & IGF1R & EGFR
trainingLabel<-list(control=list(bad=1:12,igf1r=1:12,her2=1:12,egfr=30:35),bad=13:18,igf1r=19:24,her2=25:29,egfr=36:41)
bad_igf1r_her2_egfr_genelist<-c(bad_200_genelist,igf1r_75_genelist,her2_15_genelist,egfr_25_genelist)
names(bad_igf1r_her2_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="bad_igf1r_her2_egfr", geneList=bad_igf1r_her2_egfr_genelist, trainingData = cbind(c_gfp,c_bad,c_igf1r,c_her2,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#6. HER2 & BAD & ERK & EGFR
trainingLabel<-list(control=list(bad=1:12,erk=1:12,her2=1:12,egfr=30:35),bad=13:18,erk=19:24,her2=25:29,egfr=36:41)
bad_erk_her2_egfr_genelist<-c(bad_200_genelist,erk_250_genelist,her2_15_genelist,egfr_25_genelist)
names(bad_erk_her2_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="bad_erk_her2_egfr", geneList=bad_erk_her2_egfr_genelist, trainingData = cbind(c_gfp,c_bad,c_erk,c_her2,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#7. HER2 & BAD & RAF & EGFR
trainingLabel<-list(control=list(bad=1:12,raf=1:12,her2=1:12,egfr=30:35),bad=13:18,raf=19:24,her2=25:29,egfr=36:41)
bad_raf_her2_egfr_genelist<-c(bad_200_genelist,raf_100_genelist,her2_15_genelist,egfr_25_genelist)
names(bad_raf_her2_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="bad_raf_her2_egfr", geneList=bad_raf_her2_egfr_genelist, trainingData = cbind(c_gfp,c_bad,c_raf,c_her2,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#8. HER2 & IGF1R & RAF & EGFR
trainingLabel<-list(control=list(igf1r=1:12,raf=1:12,her2=1:12,egfr=30:35),igf1r=13:18,raf=19:24,her2=25:29,egfr=36:41)
igf1r_raf_her2_egfr_genelist<-c(igf1r_75_genelist,raf_100_genelist,her2_15_genelist,egfr_25_genelist)
names(igf1r_raf_her2_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="igf1r_raf_her2_egfr", geneList=igf1r_raf_her2_egfr_genelist, trainingData = cbind(c_gfp,c_igf1r,c_raf,c_her2,train_egfr), testData = c_test, trainingLabel = trainingLabel)



#9. HER2 & IGF1R & ERK & EGFR
trainingLabel<-list(control=list(bad=1:12,raf=1:12,her2=1:12,egfr=30:35),igf1r=13:18,erk=19:24,her2=25:29,egfr=36:41)
igf1r_erk_her2_egfr_genelist<-c(igf1r_75_genelist,erk_250_genelist,her2_15_genelist,egfr_25_genelist)
names(igf1r_erk_her2_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="igf1r_erk_her2_egfr", geneList=igf1r_erk_her2_egfr_genelist, trainingData = cbind(c_gfp,c_igf1r,c_erk,c_her2,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#10. HER2 & RAF & ERK & EGFR
trainingLabel<-list(control=list(raf=1:12,erk=1:12,her2=1:12,egfr=30:35),raf=13:18,erk=19:24,her2=25:29,egfr=36:41)
raf_erk_her2_egfr_genelist<-c(raf_100_genelist,erk_250_genelist,her2_15_genelist,egfr_25_genelist)
names(raf_erk_her2_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="raf_erk_her2_egfr", geneList=raf_erk_her2_egfr_genelist, trainingData = cbind(c_gfp,c_raf,c_erk,c_her2,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#11. AKT, IGF1R, ERK & EGFR
trainingLabel<-list(control=list(akt=1:12,igf1r=1:12,erk=1:12,egfr=31:36),akt=13:18,igf1r=19:24,erk=25:30,egfr=37:42)
akt_igf1r_erk_egfr_genelist<-c(akt_75_genelist,igf1r_75_genelist,erk_250_genelist,egfr_25_genelist)
names(akt_igf1r_erk_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="akt_igf1r_erk_egfr", geneList=akt_igf1r_erk_egfr_genelist, trainingData = cbind(c_gfp,c_akt,c_igf1r,c_erk,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#12. AKT, IGF1R, RAF & EGFR
trainingLabel<-list(control=list(akt=1:12,igf1r=1:12,raf=1:12,egfr=31:36),akt=13:18,igf1r=19:24,raf=25:30,egfr=37:42)
akt_igf1r_raf_egfr_genelist<-c(akt_75_genelist,igf1r_75_genelist,raf_100_genelist,egfr_25_genelist)
names(akt_igf1r_raf_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="akt_igf1r_raf_egfr", geneList=akt_igf1r_raf_egfr_genelist, trainingData = cbind(c_gfp,c_akt,c_igf1r,c_raf,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#13. AKT, BAD, ERK & EGFR
trainingLabel<-list(control=list(akt=1:12,bad=1:12,erk=1:12,egfr=31:36),akt=13:18,bad=19:24,erk=25:30,egfr=37:42)
akt_bad_erk_egfr_genelist<-c(akt_75_genelist,bad_200_genelist,erk_250_genelist,egfr_25_genelist)
names(akt_bad_erk_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="akt_bad_erk_egfr", geneList=akt_bad_erk_egfr_genelist, trainingData = cbind(c_gfp,c_akt,c_bad,c_erk,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#14. AKT, BAD, RAF & EGFR
trainingLabel<-list(control=list(akt=1:12,bad=1:12,raf=1:12,egfr=31:36),akt=13:18,bad=19:24,raf=25:30,egfr=37:42)
akt_bad_raf_egfr_genelist<-c(akt_75_genelist,bad_200_genelist,raf_100_genelist,egfr_25_genelist)
names(akt_bad_raf_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="akt_bad_raf_egfr", geneList=akt_bad_erk_egfr_genelist, trainingData = cbind(c_gfp,c_akt,c_bad,c_raf,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#15. AKT, BAD, IGF1R & EGFR
trainingLabel<-list(control=list(akt=1:12,bad=1:12,igf1r=1:12,egfr=31:36),akt=13:18,bad=19:24,igf1r=25:30,egfr=37:42)
akt_bad_igf1r_egfr_genelist<-c(akt_75_genelist,bad_200_genelist,igf1r_75_genelist,egfr_25_genelist)
names(akt_bad_igf1r_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="akt_bad_igf1r_egfr", geneList=akt_bad_igf1r_egfr_genelist, trainingData = cbind(c_gfp,c_akt,c_bad,c_igf1r,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#16. BAD, IGF1R, ERK & EGFR
trainingLabel<-list(control=list(bad=1:12,igf1r=1:12,erk=1:12,egfr=31:36),bad=13:18,igf1r=19:24,erk=25:30,egfr=37:42)
bad_igf1r_erk_egfr_genelist<-c(bad_200_genelist,igf1r_75_genelist,erk_250_genelist,egfr_25_genelist)
names(bad_igf1r_erk_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="bad_igf1r_erk_egfr", geneList=bad_igf1r_erk_egfr_genelist, trainingData = cbind(c_gfp,c_bad,c_igf1r,c_erk,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#17. BAD, IGF1R, RAF & EGFR
trainingLabel<-list(control=list(bad=1:12,igf1r=1:12,raf=1:12,egfr=31:36),bad=13:18,igf1r=19:24,raf=25:30,egfr=37:42)
bad_igf1r_raf_egfr_genelist<-c(bad_200_genelist,igf1r_75_genelist,raf_100_genelist,egfr_25_genelist)
names(bad_igf1r_raf_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="bad_igf1r_raf_egfr", geneList=bad_igf1r_raf_egfr_genelist, trainingData = cbind(c_gfp,c_bad,c_igf1r,c_raf,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#18. BAD, RAF, ERK & EGFR
trainingLabel<-list(control=list(bad=1:12,raf=1:12,erk=1:12,egfr=31:36),bad=13:18,raf=19:24,erk=25:30,egfr=37:42)
bad_raf_erk_egfr_genelist<-c(bad_200_genelist,raf_100_genelist,erk_250_genelist,egfr_25_genelist)
names(bad_raf_erk_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="bad_raf_erk_egfr", geneList=bad_raf_erk_egfr_genelist, trainingData = cbind(c_gfp,c_bad,c_raf,c_erk,train_egfr), testData = c_test, trainingLabel = trainingLabel)

#19. RAF, IGF1R, ERK & EGFR
trainingLabel<-list(control=list(raf=1:12,igf1r=1:12,erk=1:12,egfr=31:36),raf=13:18,igf1r=19:24,erk=25:30,egfr=37:42)
raf_igf1r_erk_egfr_genelist<-c(raf_100_genelist,igf1r_75_genelist,erk_250_genelist,egfr_25_genelist)
names(raf_igf1r_erk_egfr_genelist)<-names(trainingLabel)[-1]
testSig_multi(sigProteins="raf_igf1r_erk_egfr", geneList=raf_igf1r_erk_egfr_genelist, trainingData = cbind(c_gfp,c_raf,c_igf1r,c_erk,train_egfr), testData = c_test, trainingLabel = trainingLabel)











######[,12]   [,13]   [,14]  [,15]  [,16]  [,17]   [,18]   [,19]   [,20]  [,21]   [,22]  
#[1,] "AKT"   "AKT"   "AKT"  "AKT"  "AKT"  "AKT"   "AKT"   "AKT"   "AKT"  "BAD"   "BAD"  
#[2,] "HER2"  "HER2"  "HER2" "HER2" "HER2" "IGF1R" "IGF1R" "IGF1R" "ERK"  "HER2"  "HER2" 
#[3,] "IGF1R" "IGF1R" "ERK"  "ERK"  "EGFR" "ERK"   "ERK"   "EGFR"  "EGFR" "IGF1R" "IGF1R"
#[4,] "EGFR"  "RAF"   "EGFR" "RAF"  "RAF"  "EGFR"  "RAF"   "RAF"   "RAF"  "ERK"   "EGFR" 
#[,23]   [,24]  [,25]  [,26]  [,27]   [,28]   [,29]   [,30]  [,31]   [,32]   [,33]  
#[1,] "BAD"   "BAD"  "BAD"  "BAD"  "BAD"   "BAD"   "BAD"   "BAD"  "HER2"  "HER2"  "HER2" 
#[2,] "HER2"  "HER2" "HER2" "HER2" "IGF1R" "IGF1R" "IGF1R" "ERK"  "IGF1R" "IGF1R" "IGF1R"
#[3,] "IGF1R" "ERK"  "ERK"  "EGFR" "ERK"   "ERK"   "EGFR"  "EGFR" "ERK"   "ERK"   "EGFR" 
#[4,] "RAF"   "EGFR" "RAF"  "RAF"  "EGFR"  "RAF"   "RAF"   "RAF"  "EGFR"  "RAF"   "RAF"  
######[,34]  [,35]  
#[1,] "HER2" "IGF1R"
#[2,] "ERK"  "ERK"  
#[3,] "EGFR" "EGFR" 
#[4,] "RAF"  "RAF"  
pathways_5<-combn(pathways,5)
######[,1]    [,2]    [,3]    [,4]   [,5]   [,6]   [,7]    [,8]    [,9]    [,10]  [,11]  
#[1,] "AKT"   "AKT"   "AKT"   "AKT"  "AKT"  "AKT"  "AKT"   "AKT"   "AKT"   "AKT"  "AKT"  
#[2,] "BAD"   "BAD"   "BAD"   "BAD"  "BAD"  "BAD"  "BAD"   "BAD"   "BAD"   "BAD"  "HER2" 
#[3,] "HER2"  "HER2"  "HER2"  "HER2" "HER2" "HER2" "IGF1R" "IGF1R" "IGF1R" "ERK"  "IGF1R"
#[4,] "IGF1R" "IGF1R" "IGF1R" "ERK"  "ERK"  "EGFR" "ERK"   "ERK"   "EGFR"  "EGFR" "ERK"  
#[5,] "ERK"   "EGFR"  "RAF"   "EGFR" "RAF"  "RAF"  "EGFR"  "RAF"   "RAF"   "RAF"  "EGFR" 
#######[,12]   [,13]   [,14]  [,15]   [,16]   [,17]   [,18]   [,19]  [,20]   [,21]  
#[1,] "AKT"   "AKT"   "AKT"  "AKT"   "BAD"   "BAD"   "BAD"   "BAD"  "BAD"   "HER2" 
#[2,] "HER2"  "HER2"  "HER2" "IGF1R" "HER2"  "HER2"  "HER2"  "HER2" "IGF1R" "IGF1R"
#[3,] "IGF1R" "IGF1R" "ERK"  "ERK"   "IGF1R" "IGF1R" "IGF1R" "ERK"  "ERK"   "ERK"  
#[4,] "ERK"   "EGFR"  "EGFR" "EGFR"  "ERK"   "ERK"   "EGFR"  "EGFR" "EGFR"  "EGFR" 
#[5,] "RAF"   "RAF"   "RAF"  "RAF"   "EGFR"  "RAF"   "RAF"   "RAF"  "RAF"   "RAF" 
pathways_6<-combn(pathways,6)
######[,1]    [,2]    [,3]    [,4]   [,5]    [,6]    [,7]   
#[1,] "AKT"   "AKT"   "AKT"   "AKT"  "AKT"   "AKT"   "BAD"  
#[2,] "BAD"   "BAD"   "BAD"   "BAD"  "BAD"   "HER2"  "HER2" 
#[3,] "HER2"  "HER2"  "HER2"  "HER2" "IGF1R" "IGF1R" "IGF1R"
#[4,] "IGF1R" "IGF1R" "IGF1R" "ERK"  "ERK"   "ERK"   "ERK"  
#[5,] "ERK"   "ERK"   "EGFR"  "EGFR" "EGFR"  "EGFR"  "EGFR" 
#[6,] "EGFR"  "RAF"   "RAF"   "RAF"  "RAF"   "RAF"   "RAF"  
print(paste("Total number of pathway combinations to consider:",7+21+35+35+21+7+1,sep=" "))
dir.create("~/Dropbox/bild_signatures/gene_list_ICBP_7")
basedir<-"~/Dropbox/bild_signatures/gene_list_ICBP_7"
###########
#5-2

#5-4

#5-7

#5-9

#5-10

#5-11

#5-13

#5-14

#5-15

#5-16

#5-18

#5-19

#5-20

#5-21



###########
#6-1
trainingLabel<-list(control=list(egfr=1:6,akt=13:24,bad=13:24,her2=13:24,igf1r=13:24,erk=13:24),egfr=7:12,akt=25:30,bad=31:36,her2=37:41,igf1r=42:47,erk=48:53)
dir.create(paste(basedir,"egfr_akt_bad_her2_igf1r_erk","6_1",sep='/'))
assign_easy_multi(trainingData = cbind(train_egfr,c_gfp,c_akt,c_bad,c_her2,c_igf1r,c_erk),testData = c_test,trainingLabel1 = trainingLabel,geneList=c(egfr_genelist,akt_genelist,bad_genelist,her2_genelist,igf1r_genelist,erk_genelist),out_dir_base = paste(basedir,"egfr_akt_bad_her2_igf1r_erk","6_1",sep='/'))
#6-2
trainingLabel<-list(control=list(akt=1:12,bad=1:12,igf1r=1:12,her2=1:12,raf=1:12,erk=1:12),akt=13:18, bad=19:24,igf1r=25:30,her2=31:35, raf=36:41,erk=42:47)
sub_dir=paste(basedir,"akt_bad_igf1r_her2_raf_erk","6_1",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = cbind(c_gfp,c_akt,c_bad,c_igf1r,c_her2,c_raf,c_erk),test=c_test,trainingLabel1 = trainingLabel,geneList=c(akt_genelist,bad_genelist,igf1r_genelist,her2_genelist,raf_genelist,erk_genelist),out_dir_base = sub_dir)
#6-3
trainingLabel<-list(control=list(egfr=1:6,akt=13:24,bad=13:24,her2=13:24,igf1r=13:24,raf=13:24),egfr=7:12,akt=25:30,bad=31:36,her2=37:41,igf1r=42:47,raf=48:53)
sub_dir<-paste(basedir,"egfr_akt_bad_her2_igf1r_raf",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = cbind(train_egfr,c_gfp,c_akt,c_bad,c_her2,c_igf1r,c_raf),testData = c_test,trainingLabel1 = trainingLabel,geneList=c(egfr_genelist,akt_genelist,bad_genelist,her2_genelist,igf1r_genelist,raf_genelist),out_dir_base = paste(basedir,"egfr_akt_bad_her2_igf1r_raf",sep='/'))
#6-4
trainingLabel<-list(control=list(egfr=1:6,akt=13:24,bad=13:24,her2=13:24,raf=13:24,erk=13:24),egfr=7:12,akt=25:30,bad=31:36,her2=37:41,raf=42:47,erk=48:53)
sub_dir<-paste(basedir,"egfr_akt_bad_her2_raf_erk",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = cbind(train_egfr,c_gfp,c_akt,c_bad,c_her2,c_raf,c_erk),testData = c_test,trainingLabel1 = trainingLabel,geneList=c(egfr_genelist,akt_genelist,bad_genelist,her2_genelist,raf_genelist,erk_genelist),out_dir_base = sub_dir)

#6-5
trainingLabel<-list(control=list(egfr=1:6,akt=13:24,bad=13:24,igf1r=13:24,raf=13:24,erk=13:24),egfr=7:12,akt=25:30,bad=31:36,igf1r=37:42,raf=43:48,erk=49:54)
sub_dir<-paste(basedir,"egfr_akt_bad_igf1r_raf_erk",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = cbind(train_egfr,c_gfp,c_akt,c_bad,c_igf1r,c_raf,c_erk),testData = c_test,trainingLabel1 = trainingLabel,geneList=c(egfr_genelist,akt_genelist,bad_genelist,igf1r_genelist,raf_genelist,erk_genelist),out_dir_base = sub_dir)

#6-6
trainingLabel<-list(control=list(egfr=1:6,akt=13:24,her2=13:24,igf1r=13:24,raf=13:24,erk=13:24),egfr=7:12,akt=25:30,her2=31:35,igf1r=36:41,raf=42:47,erk=48:53)
sub_dir<-paste(basedir,"egfr_akt_her2_igfr_raf_erk",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = cbind(train_egfr,c_gfp,c_akt,c_her2,c_igf1r,c_raf,c_erk),testData = c_test,trainingLabel1 = trainingLabel,geneList=c(egfr_genelist,akt_genelist,her2_genelist,igf1r_genelist,raf_genelist,erk_genelist),out_dir_base = sub_dir)

#6-7
trainingLabel<-list(control=list(egfr=1:6,bad=13:24,her2=13:24,igf1r=13:24,raf=13:24,erk=13:24),egfr=7:12,bad=25:30,her2=31:35,igf1r=36:41,raf=42:47,erk=48:53)
sub_dir<-paste(basedir,"egfr_bad_her2_igfr_raf_erk",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = cbind(train_egfr,c_gfp,c_bad,c_her2,c_igf1r,c_raf,c_erk),testData = c_test,trainingLabel1 = trainingLabel,geneList=c(egfr_genelist,bad_genelist,her2_genelist,igf1r_genelist,raf_genelist,erk_genelist),out_dir_base = sub_dir)

#7-all

trainingLabel<-list(control=list(egfr=1:6,akt=13:24,bad=13:24,her2=13:24,igf1r=13:24,raf=13:24,erk=13:24),egfr=7:12,akt=25:30,bad=31:36,her2=37:41,igf1r=42:47,raf=48:53,erk=54:59)

#assign.wrapper(trainingData=combat_expr1[,1:59], testData=c_test,trainingLabel=trainingLabel, testLabel=NULL,geneList=c(egfr_genelist,akt_genelist,bad_genelist,her2_genelist,igf1r_genelist,raf_genelist,erk_genelist), n_sigGene=NULL, adaptive_B=TRUE,adaptive_S=FALSE, mixture_beta=TRUE,outputDir= ".", iter=2000, burn_in=1000)
assign_easy_multi(trainingData = combat_expr1[,1:59],testData = c_test,trainingLabel1 = trainingLabel,geneList=c(egfr_genelist,akt_genelist,bad_genelist,her2_genelist,igf1r_genelist,raf_genelist,erk_genelist),out_dir_base = paste(basedir,"egfr_akt_bad_her2_igf1r_erk_raf",sep='/'))
