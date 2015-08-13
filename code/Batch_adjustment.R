library(devtools)
install_github("wevanjohnson/ASSIGN")
library(ASSIGN)
install_github("wevanjohnson/sva-devel")
library(sva)


source("~/Dropbox/bild_signatures/bild_signatures/Rmarkdowns_scripts//Key_ASSIGN_functions.Rmd")
setwd("~/Dropbox/Datasets/")
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
gfp_kras<-read.table("36_GFP_KRAS_TPMlog2.txt", sep='\t', header=1, row.names=1)
#head(gfp_kras)
gfp_egfr_kras_multi_f<-merge_drop(gfp_egfr_multi_f,gfp_kras)
dim(gfp_egfr_kras_multi_f)
icbp<-as.matrix(read.table("~/Dropbox/bild_signatures/Datasets/icbp_Rsubread_tpmlog.txt", sep='\t', stringsAsFactors=FALSE, header=1, row.names=1))
all_icbp_f<-merge_drop(gfp_egfr_kras_multi_f,icbp,by=0)
colnames(all_icbp_f)
test<-all_icbp_f[,96:150]
sub<-c(6,6,12,6,6,5,6,6,6,9,9,9,9,ncol(test))
pcaplot(all_icbp_f,sub)
bat1<-as.matrix(cbind(colnames(all_icbp_f),c(rep(1,ncol(control_egfr_l)),rep(2,ncol(expr_all_f)),rep(3,ncol(gfp_kras)),rep(4,ncol(test))),c(rep(1,6),rep(2,6),rep(1,12),rep(3,6),rep(4,6),rep(5,5),rep(6,6),rep(7,6),rep(8,6),rep(1,9),rep(9,9),rep(10,9),rep(11,9),rep(1,55))))
#bat1
###one step normalization method
combat_expr1<-ComBat(dat=all_icbp_f,batch=bat1[,2], mod=NULL)#, numCovs=NULL)

###one step with GFP and ICBP with same covariate and mean only ComBat
#conditions<-c(rep("Control",6),rep("EGFR",6),rep("Control",12),rep("AKT",6),rep("BAD",6),rep("HER2",5),rep("IGF1R",6),rep("RAF",6),rep("ERK",6),rep("Control",9),rep("KRASWT",9),rep("KRASQH",9),rep("KRASGV",9),rep("ICBP",55))
mod <- model.matrix(~as.factor(bat1[,3]))
combat_expr1_mean<-ComBat(dat=all_icbp_f,batch=bat1[,2], mod=mod,mean.only = T)#, numCovs=NULL)
pcaplot(combat_expr1_mean,sub)


