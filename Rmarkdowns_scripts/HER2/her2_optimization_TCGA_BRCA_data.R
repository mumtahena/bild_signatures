set.seed(1234)
library(sva)
library(ASSIGN)
setwd("~/Dropbox/EVAN/9_17/")
feature<-as.matrix(read.table("hmec_gfp_akt_bad_gfp48_p110_rpkmlog.txt", sep='\t', stringsAsFactors=FALSE, header=1, row.names=1))
head(feature)
tcga<-as.matrix(read.table("TCGA_PANCAN20_Rsubread_BRCA_RPKMlog.txt", sep='\t', stringsAsFactors=FALSE, header=1, row.names=1,check.names=F))
head(tcga)
gfp_her2<-subset(feature,select=c(GFP.1:GFP.12,HER2.1:HER2.6))
head(gfp_her2)
nrow(gfp_her2)
gfp_her2_1 <- gfp_her2[apply(gfp_her2[,1:17]==0,1,mean) < 0.8,]###if not removed over 83%; error occurs in batch adjustment
nrow(gfp_her2_1)
AVG_IN <- apply(gfp_her2_1[,1:17], 1, mean)
VAR_IN <- apply(gfp_her2_1[,1:17], 1, var)
threshold_E <- sort(AVG_IN)[round(length(AVG_IN)*0.1)]
threshold_V <- sort(VAR_IN)[round(length(VAR_IN)*0.1)]
filtered_feature <- gfp_her2_1[which((AVG_IN > threshold_E)&(VAR_IN > threshold_V)),]
nrow(filtered_feature)
nrow(gfp_her2)
nrow(gfp_her2_1)
gfp_her2_1_tcga<-merge(filtered_feature,tcga,by=0)
nrow(gfp_her2_1_tcga)
rownames(gfp_her2_1_tcga)<-gfp_her2_1_tcga$Row.names
gfp_her2_1_tcga<-gfp_her2_1_tcga[,2:ncol(gfp_her2_1_tcga)]

#running ASSIGN w/o any batch normalization
####################
trainingLabel1 <- list(control = list(control1=1:12), her2=13:17) 

####for TCGA BRCA samples####
ncol(gfp_her2_1_tcga)
train_her2_2<-gfp_her2_1_tcga[,1:17]
test_nocombat2<-gfp_her2_1_tcga[,18:1099]
setwd("~/Dropbox/EVAN/9_17/Results/TCGA/")
dir.create("her2_single_100_nocombat_tcga_f/")#, 
set.seed(1234)
assign.wrapper(trainingData=train_her2_2, testData=test_nocombat2, trainingLabel=trainingLabel1, testLabel=NULL, geneList=NULL, n_sigGene=100, adaptive_B=F, adaptive_S=F, mixture_beta=F, outputDir="her2_single_100_nocombat_tcga_f/", theta0=0.05, theta1=0.9, iter=2000, burn_in=1000)
dir.create("her2_single_Adap_100_nocombat_tcga_f/")#,
set.seed(1234)
assign.wrapper(trainingData=train_her2_2, testData=test_nocombat2, trainingLabel=trainingLabel1, testLabel=NULL, geneList=NULL, n_sigGene=100, adaptive_B=F, adaptive_S=F, mixture_beta=F, outputDir="her2_single_Adap_100_nocombat_tcga_f/", theta0=0.05, theta1=0.9, iter=2000, burn_in=1000)
##########combat with cov

head(gfp_her2_1_tcga)
pca_combat <- prcomp(t(gfp_her2_1_tcga), center=T,scale=T)
plot(pca_combat)
plot(pca_combat$x[,1],pca_combat$x[,2])
ncol(pca_combat)
points(pca_combat$x[1:12],pca_combat$x[1:12,2],col=2)#gfp in red
points(pca_combat$x[13:17],pca_combat$x[13:17,2],col=4)#her2 in blue
points(pca_combat$x[18:1099,1],pca_combat$x[18:1099,2],col=3)#tcga greens
#with covariates
bat1<-as.matrix(cbind(c(colnames(gfp_her2_1),colnames(tcga)),c(rep(1,length(colnames(gfp_her2_1))),rep(2,length(colnames(tcga)))),c(rep(1,12),rep(2,5),rep(1,ncol(tcga)))))
head(bat1)
bat1[,2]
mod <- model.matrix(~as.factor(bat1[,3]))
##running combat with covariates
combat_expr1<-ComBat(dat=gfp_her2_1_tcga, batch=bat1[,2], mod=mod, numCovs=NULL)
pca_combat <- prcomp(t(combat_expr1), center=T,scale=T)
plot(pca_combat)
plot(pca_combat$x[,1],pca_combat$x[,2])
points(pca_combat$x[1:12],pca_combat$x[1:12,2],col=2)#gfp in red
points(pca_combat$x[13:17],pca_combat$x[13:17,2],col=4)#her2 in blue
points(pca_combat$x[18:1099,1],pca_combat$x[18:1099,2],col=3)#tcga greens
ncol(combat_expr1)
write.table(combat_expr1[,1:12],file="post_combat_gfp_cov_tcga_f.txt",col.names=NA,sep='\t',quote=F)
write.table(combat_expr1[,13:17],file="post_combat_her2_cov_tcga_f.txt",col.names=NA,sep='\t',quote=F)
write.table(combat_expr1[,18:1099],file="post_combat_cov_tcga_f.txt",col.names=NA,sep='\t',quote=F)

train_her2<-combat_expr1[,1:17]
test<-combat_expr1[,18:1099]
head(train_her2)
nrow(train_her2)

dir.create("her2_single_100_adap_cov_tcga_f/")#, 
set.seed(1234)
assign.wrapper(trainingData=train_her2, testData=test, trainingLabel=trainingLabel1, testLabel=NULL, geneList=NULL, n_sigGene=100, adaptive_B=T, adaptive_S=F, mixture_beta=F, outputDir="her2_single_100_adap_cov_tcga_f/", theta0=0.05, theta1=0.9, iter=2000, burn_in=1000)

dir.create("her2_single_100_nonadap_cov_tcga_f/")#, 
set.seed(1234)
assign.wrapper(trainingData=train_her2, testData=test, trainingLabel=trainingLabel1, testLabel=NULL, geneList=NULL, n_sigGene=100, adaptive_B=F, adaptive_S=F, mixture_beta=F, outputDir="her2_single_100_nonadap_cov_tcga_f/", theta0=0.05, theta1=0.9, iter=2000, burn_in=1000)


############combat without cov

bat2<-as.matrix(cbind(c(colnames(gfp_her2_1),colnames(tcga)),c(rep(1,length(colnames(gfp_her2_1))),rep(2,length(colnames(tcga))))))
head(bat2)
bat2[,2]
##running combat without covariates
combat_expr1<-ComBat(dat=gfp_her2_1_tcga, batch=bat2[,2], mod=NULL, numCovs=NULL)
pca_combat <- prcomp(t(combat_expr1), center=T,scale=T)
plot(pca_combat)
plot(pca_combat$x[,1],pca_combat$x[,2])
points(pca_combat$x[1:12],pca_combat$x[1:12,2],col=2)#gfp in red
points(pca_combat$x[13:17],pca_combat$x[13:17,2],col=4)#her2 in blue
points(pca_combat$x[18:1099,1],pca_combat$x[18:1099,2],col=3)#tcga greens
ncol(combat_expr1)
write.table(combat_expr1[,1:12],file="post_combat_gfp_wo_cov_tcga_f.txt",col.names=NA,sep='\t',quote=F)
write.table(combat_expr1[,13:17],file="post_combat_her2_wo_cov_tcga_f.txt",col.names=NA,sep='\t',quote=F)
write.table(combat_expr1[,18:1099],file="post_combat_wo_cov_tcga_f.txt",col.names=NA,sep='\t',quote=F)

train_her2<-combat_expr1[,1:17]
test<-combat_expr1[,18:1099]
head(train_her2)
nrow(train_her2)

dir.create("her2_single_100_adap_wo_cov_tcga_f/")#, 
set.seed(1234)
assign.wrapper(trainingData=train_her2, testData=test, trainingLabel=trainingLabel1, testLabel=NULL, geneList=NULL, n_sigGene=100, adaptive_B=T, adaptive_S=F, mixture_beta=F, outputDir="her2_single_100_adap_wo_cov_tcga_f/", theta0=0.05, theta1=0.9, iter=2000, burn_in=1000)

dir.create("her2_single_100_nonadap_wo_cov_tcga_f/")#, 
set.seed(1234)
assign.wrapper(trainingData=train_her2, testData=test, trainingLabel=trainingLabel1, testLabel=NULL, geneList=NULL, n_sigGene=100, adaptive_B=F, adaptive_S=F, mixture_beta=F, outputDir="her2_single_100_nonadap_wo_cov_tcga_f/", theta0=0.05, theta1=0.9, iter=2000, burn_in=1000)

############ASSIGN after dwd-bild###########
dwd_bild1<-read.table("~/Dropbox/EVAN/9_17/Results/TCGA/her2_opt_f/attic/dataset.dwd_bild.gct", sep='\t', header=1, check.names =F,skip = 2)
rownames(dwd_bild1)<-dwd_bild1$NAME
dwd_bild1<-subset(dwd_bild1,select=-c(NAME,DESCRIPTION)
pca_combat <- prcomp(t(dwd_bild1), center=T,scale=T)
plot(pca_combat)
plot(pca_combat$x[,1],pca_combat$x[,2])
points(pca_combat$x[1:12],pca_combat$x[1:12,2],col=2)#gfp in red
points(pca_combat$x[13:17],pca_combat$x[13:17,2],col=4)# her2 in blue
points(pca_combat$x[18:1099,1],pca_combat$x[18:1099,2],col=3)#icbp greens

trainingLabel1 <- list(control = list(control1=1:12), her2=13:17) 
trainingLabel1
train_her2<-dwd_bild1[,1:17]
test<-dwd_bild1[,18:1099]
head(train_her2)
nrow(train_her2)
dir.create("her2_single_100_dwd_adap_tcga_f/")#, 
set.seed(1234)
assign.wrapper(trainingData=train_her2, testData=test, trainingLabel=trainingLabel1, testLabel=NULL, geneList=NULL, n_sigGene=100, adaptive_B=T, adaptive_S=F, mixture_beta=F, outputDir="her2_single_100_dwd_adap_tcga_f/", theta0=0.05, theta1=0.9, iter=2000, burn_in=1000)

dir.create("her2_single_100_dwd_non_adap_tcga_f/")#, 
set.seed(1234)
assign.wrapper(trainingData=train_her2, testData=test, trainingLabel=trainingLabel1, testLabel=NULL, geneList=NULL, n_sigGene=100, adaptive_B=F, adaptive_S=F, mixture_beta=F, outputDir="her2_single_100_dwd_non_adap_tcga_f/", theta0=0.05, theta1=0.9, iter=2000, burn_in=1000)

##########################
setwd("~/Dropbox/EVAN/9_17/")
cnv<-read.table("cnv_erbb2_tumor", sep='\t',row.names=1, header=1, check.names = F)
##converting to shorter barcode ids
a<-NULL
tmp<-NULL
a<-strsplit(rownames(cnv),'-')
for(i in 1:length(a)){tmp[[i]]<-paste(a[[i]][1:3],collapse='-')}
#head(tmp)
#tmp
for(i in 1:length(a)){cnv[i,"truncated"]<-tmp[[i]][1]}
head(cnv)

results<-read.table("10_6.txt", header=1, row.names=1)
head(results)

##converting to shorter barcode ids
a<-NULL
tmp<-NULL
a<-strsplit(rownames(results),'-')
for(i in 1:length(a)){tmp[[i]]<-paste(a[[i]][1:3],collapse='-')}
head(tmp)
#tmp
for(i in 1:length(a)){results[i,"truncated"]<-tmp[[i]][1]}
head(results)
final_mat<-merge(results,cnv,by.x='truncated',by.y='truncated')
head(final_mat)
sum(duplicated(final_mat$truncated))

###checking ASSIGN outputs######
cor(final_mat$A_no_bat_NonAdap,final_mat$ERBB2)
cor(final_mat$A_no_bat_Adap,final_mat$ERBB2)
cor(final_mat$A_combat_cov_NonAdap,final_mat$ERBB2)
cor(final_mat$A_combat_cov_adap,final_mat$ERBB2)
cor(final_mat$A_combat_wo_cov_NonAdap,final_mat$ERBB2)
cor(final_mat$A_combat_wo_cov_adap,final_mat$ERBB2)
cor(final_mat$A_combat_wo_cov_NonAdap_nofilter,final_mat$ERBB2)
cor(final_mat$A_combat_wo_cov_adap_no_filter,final_mat$ERBB2)
cor(final_mat$A_dwd_bild_NonAdaptive_f,final_mat$ERBB2)
cor(final_mat$A_dwd_bild_Adaptive_f,final_mat$ERBB2)

#####checking ASSIGN No Combat output###
cor(final_mat$B_dwd_bild_f,final_mat$ERBB2)
cor(final_mat$B_dwd_bild,final_mat$ERBB2)
cor(final_mat$B_combat_wo_cov_f,final_mat$ERBB2)
cor(final_mat$B_combat_cov_f,final_mat$ERBB2)
#####checking ASSIGN Combat with covariate output###
cor(final_mat$cov_ComBat_NA,final_mat$ERBB2)
cor(final_mat$cov_ComBat_A,final_mat$ERBB2)
#####checking ASSIGN simple Combat output###
cor(final_mat$ComBat_NA,final_mat$ERBB2)
cor(final_mat$ComBat_A,final_mat$ERBB2)
#####checking dwd_bild ASSIGN output###
cor(final_mat$DWD_ASSIGN_NA,final_mat$ERBB2)
cor(final_mat$DWD_ASSIGN_A,final_mat$ERBB2)


