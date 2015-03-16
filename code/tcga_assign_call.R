setwd("~/Dropbox/bild_signatures/Datasets/")
expr<-as.matrix(read.table("GFP18_AKT_BAD_HER2_IGF1R_RAF_ERK.tpmlog",sep='\t',row.names=1,header=1))
control<-subset(expr, select=GFP.1:GFP.12)
her2<-subset(expr, select=HER2.1:HER2.6)
akt<-subset(expr,select=AKT.1:AKT.6)
bad<-subset(expr,select=BAD.1:BAD.6)
igf1r<-subset(expr,select=IGF1R.1:IGF1R.6)
erk<-subset(expr,select=ERK.1:ERK.6)
expr_all<-cbind(control,akt,bad,her2,igf1r,erk)
tcga<-as.matrix(read.table("~/Dropbox/Datasets/TCGA20_brca_1_23.txt", sep='\t', stringsAsFactors=T, header=1, row.names=1))
expr_all_f <-expr_all[apply(expr_all[,1:41]==0,1,mean) < 0.85,]
dim(expr_all_f)
expr_all_tcga_f<-merge_drop(expr_all_f,tcga,by=0)
dim(expr_all_tcga_f)
sub<-c(12,6,6,5,6,6,length(colnames(tcga)))
pdf(file='~/Dropbox/bild_signatures/tcga_hmec_pca_plot_3_14_15.pdf')
pcaplot(expr_all_tcga_f,sub)
bat1<-as.matrix(cbind(c(colnames(expr_all_f),colnames(tcga)),c(rep(1,length(colnames(expr_all_f))),rep(2,length(colnames(tcga))))))
#bat1
combat_expr1<-ComBat(dat=expr_all_tcga_f, batch=bat1[,2], mod=NULL, numCovs=NULL)
pcaplot(combat_expr1,sub)
dev.off()
combat_tcga<-combat_expr1;#as.matrix(read.table("~/Dropbox/Datasets/TCGA20_brca_hmec_combat.txt", sep='\t', stringsAsFactors=T, header=1, row.names=1))
c_gfp<-subset(combat_tcga, select=GFP.1:GFP.12)
c_akt<-subset(combat_tcga, select=AKT.1:AKT.6)
c_bad<-subset(combat_tcga, select=BAD.1:BAD.6)
c_her2<-subset(combat_tcga, select=HER2.1:HER2.6)
c_igf1r<-subset(combat_tcga, select=IGF1R.1:IGF1R.6)
c_erk<-subset(combat_tcga, select=ERK.1:ERK.6)
c_test<-combat_tcga[,42:ncol(combat_tcga)]
basedir="~/Dropbox/bild_signatures/tcga_15_mar_all"
dir.create( basedir)

#############trying one pathways at a time in multipathway#############
#1. AKT
trainingLabela<-list(control=list(akt=1:12),akt=13:18)
sub_dir<-paste(basedir,"akt",sep='/')
dir.create( sub_dir)
assign_easy_multi(trainingData = cbind(c_gfp,c_akt),test=c_test,trainingLabel1 = trainingLabela,g=150,out_dir_base = sub_dir,single = 1)

#2. BAD
trainingLabelb<-list(control=list(bad=1:12),bad=13:18)
sub_dir<-paste(basedir,"bad",sep='/')
dir.create( sub_dir)
assign_easy_multi(trainingData = cbind(c_gfp,c_bad),test=c_test,trainingLabel1 = trainingLabelb,g=150,out_dir_base = sub_dir,single = 1)

#3. HER2
trainingLabelh<-list(control=list(her2=1:12),her2=13:17)
sub_dir<-paste(basedir,"her2",sep='/')
dir.create( sub_dir)
assign_easy_multi(trainingData = cbind(c_gfp,c_her2),test=c_test,trainingLabel1 = trainingLabelh,g=15,out_dir_base = sub_dir,single = 1)

#4. IGF1R
trainingLabeli<-list(control=list(igf1r=1:12),igf1r=13:18)
sub_dir<-paste(basedir,"igf1r",sep='/')
dir.create( sub_dir)
assign_easy_multi(trainingData = cbind(c_gfp,c_igf1r),test=c_test,trainingLabel1 = trainingLabeli,g=100,out_dir_base = sub_dir,single = 1)

#5. ERK
trainingLabele<-list(control=list(erk=1:12),erk=13:18)
sub_dir<-paste(basedir,"erk",sep='/')
dir.create( sub_dir)
assign_easy_multi(trainingData = cbind(c_gfp,c_erk),test=c_test,trainingLabel1 = trainingLabele,g=100,out_dir_base = sub_dir,single = 1)

#############trying two pathways at a time in multipathway#############
#1. HER2 & AKT
trainha<-cbind(c_gfp,c_akt,c_her2)
trainingLabelha<-list(control=list(akt=1:12,her2=1:12),akt=13:18,her2=19:23)
sub_dir=paste(basedir,"her2_akt",sep='/')
dir.create( sub_dir)
assign_easy_multi(trainingData = trainha,test=c_test,trainingLabel1 = trainingLabelha,g=c(150,15),out_dir_base = sub_dir)

#2. HER2 & BAD
trainhb<-cbind(c_gfp,c_bad,c_her2)
trainingLabelhb<-list(control=list(bad=1:12,her2=1:12),bad=13:18,her2=19:23)
sub_dir=paste(basedir,"her2_bad",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = trainhb,test=c_test,trainingLabel1 = trainingLabelhb,g=c(150,15),out_dir_base = sub_dir)

#3. HER2 & IGF1R
trainhi<-cbind(c_gfp,c_igf1r,c_her2)
trainingLabelhi<-list(control=list(igf1r=1:12,her2=1:12),igf1r=13:18,her2=19:23)
sub_dir=paste(basedir,"her2_igf1r",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = trainhi,test=c_test,trainingLabel1 = trainingLabelhi,g=c(100,15),out_dir_base = sub_dir)

#4. AKT & BAD
trainab<-cbind(c_gfp,c_akt,c_bad)
trainingLabelab<-list(control=list(akt=1:12,bad=1:12),akt=13:18,bad=19:24)
sub_dir=paste(basedir,"akt_bad",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = trainab,test=c_test,trainingLabel1 = trainingLabelab,g=c(150,150),out_dir_base = sub_dir)

#5. AKT & IGF1R
trainai<-cbind(c_gfp,c_akt,c_igf1r)

trainingLabelai<-list(control=list(akt=1:12,igf1r=1:12),akt=13:18,igf1r=19:24)
sub_dir=paste(basedir,"akt_igf1r",sep='/')
dir.create( sub_dir)
assign_easy_multi(trainingData = trainai,test=c_test,trainingLabel1 = trainingLabelai,g=c(150,100),out_dir_base = sub_dir)

#6. BAD & IGF1R
trainbi<-cbind(c_gfp,c_bad,c_igf1r)
trainingLabelbi<-list(control=list(bad=1:12,igf1r=1:12),bad=13:18,igf1r=19:24)
sub_dir=paste(basedir,"bad_igf1r",sep='/')
dir.create( sub_dir)
assign_easy_multi(trainingData = trainbi,test=c_test,trainingLabel1 = trainingLabelbi,g=c(150,100),out_dir_base = sub_dir)

#7. ERK & IGF1R
trainei<-cbind(c_gfp,c_erk,c_igf1r)
trainingLabelei<-list(control=list(erk=1:12,igf1r=1:12),erk=13:18,igf1r=19:24)
sub_dir=paste(basedir,"erk_igf1r",sep='/')
dir.create( sub_dir)
assign_easy_multi(trainingData = trainei,test=c_test,trainingLabel1 = trainingLabelei,g=c(100,100),out_dir_base = sub_dir)

#8. ERK & AKT
trainea<-cbind(c_gfp,c_erk,c_akt)
trainingLabelea<-list(control=list(erk=1:12,akt=1:12),erk=13:18,akt=19:24)
sub_dir=paste(basedir,"erk_akt",sep='/')
dir.create( sub_dir)
assign_easy_multi(trainingData = trainea,test=c_test,trainingLabel1 = trainingLabelea,g=c(100,150),out_dir_base = sub_dir)

#9. ERK & BAD
traineb<-cbind(c_gfp,c_erk,c_bad)
trainingLabeleb<-list(control=list(erk=1:12,bad=1:12),erk=13:18,bad=19:24)
sub_dir=paste(basedir,"erk_bad",sep='/')
dir.create( sub_dir)
assign_easy_multi(trainingData = traineb,test=c_test,trainingLabel1 = trainingLabeleb,g=c(100,150),out_dir_base = sub_dir)

#10. ERK & HER2
traineh<-cbind(c_gfp,c_erk,c_her2)
trainingLabeleh<-list(control=list(erk=1:12,her2=1:12),erk=13:18,her2=19:23)
sub_dir=paste(basedir,"erk_her2",sep='/')
dir.create( sub_dir)
assign_easy_multi(trainingData = traineh,test=c_test,trainingLabel1 = trainingLabeleh,g=c(100,15),out_dir_base = sub_dir)

#############trying three pathways at a time in multipathway#############
#1. HER2, AKT & BAD
trainhab<-cbind(c_gfp,c_akt,c_bad,c_her2)
trainingLabelhab<-list(control=list(akt=1:12,bad=1:12,her2=1:12),akt=13:18,bad=19:24,her2=25:29)
sub_dir=paste(basedir,"her2_akt_bad",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = trainhab,test=c_test,trainingLabel1 = trainingLabelhab,g=c(150,150,15),out_dir_base = sub_dir)

#2. HER2, BAD & IGF1R
trainhbi<-cbind(c_gfp,c_igf1r,c_bad,c_her2)
trainingLabelhbi<-list(control=list(igf1r=1:12,bad=1:12,her2=1:12),igf1r=13:18,bad=19:24,her2=25:29)
sub_dir=paste(basedir,"her2_bad_igf1r",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = trainhbi,test=c_test,trainingLabel1 = trainingLabelhbi,g=c(100,150,15),out_dir_base = sub_dir)

#3. AKT, BAD & IGF1R
trainabi<-cbind(c_gfp,c_akt,c_bad,c_igf1r)
trainingLabelabi<-list(control=list(akt=1:12,bad=1:12,igf1r=1:12),akt=13:18,bad=19:24,igf1r=25:30)
sub_dir=paste(basedir,"akt_bad_igf1r",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = trainabi,test=c_test,trainingLabel1 = trainingLabelabi,g=c(150,150,100),out_dir_base = sub_dir)

#4. AKT, BAD & ERK
trainabe<-cbind(c_gfp,c_akt,c_bad,c_erk)
trainingLabelabe<-list(control=list(akt=1:12,bad=1:12,erk=1:12),akt=13:18,bad=19:24,erk=25:30)
sub_dir=paste(basedir,"akt_bad_erk",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = trainabe,test=c_test,trainingLabel1 = trainingLabelabe,g=c(150,150,100),out_dir_base = sub_dir)

#5. AKT, HER2 & IGF1R
trainahi<-cbind(c_gfp,c_akt,c_her2,c_igf1r)
trainingLabelahi<-list(control=list(akt=1:12,her2=1:12,igf1r=1:12),akt=13:18,her2=19:23,igf1r=24:29)
sub_dir=paste(basedir,"akt_her2_igf1r",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = trainahi,test=c_test,trainingLabel1 = trainingLabelahi,g=c(150,15,100),out_dir_base = sub_dir)

#6. AKT, HER2 & ERK
trainahe<-cbind(c_gfp,c_akt,c_her2,c_erk)
trainingLabelahe<-list(control=list(akt=1:12,her2=1:12,igf1r=1:12),akt=13:18,her2=19:23,erk=24:29)
sub_dir=paste(basedir,"akt_her2_erk",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = trainahe,test=c_test,trainingLabel1 = trainingLabelahe,g=c(150,15,100),out_dir_base = sub_dir)

#7. AKT, IGF1R & ERK
trainaie<-cbind(c_gfp,c_akt,c_igf1r,c_erk)
trainingLabelaie<-list(control=list(akt=1:12,igf1r=1:12,erk=1:12),akt=13:18,igf1r=19:24,erk=25:30)
sub_dir=paste(basedir,"akt_igf1r_erk",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = trainaie,test=c_test,trainingLabel1 = trainingLabelaie,g=c(150,100,100),out_dir_base = sub_dir)

#8. BAD, IGF1R & ERK
trainbie<-cbind(c_gfp,c_bad,c_igf1r,c_erk)
trainingLabelbie<-list(control=list(bad=1:12,igf1r=1:12,erk=1:12),bad=13:18,igf1r=19:24,erk=25:30)
sub_dir=paste(basedir,"bad_igf1r_erk",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = trainbie,test=c_test,trainingLabel1 = trainingLabelbie,g=c(150,100,100),out_dir_base = sub_dir)

#9. BAD, HER2 & ERK
trainbhe<-cbind(c_gfp,c_bad,c_her2,c_erk)
trainingLabelbhe<-list(control=list(bad=1:12,her2=1:12,erk=1:12),bad=13:18,her2=19:23,erk=24:29)
sub_dir=paste(basedir,"bad_her2_erk",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = trainbhe,test=c_test,trainingLabel1 = trainingLabelbhe,g=c(150,15,100),out_dir_base = sub_dir)

#10. IGF1R, HER2 & ERK
trainihe<-cbind(c_gfp,c_igf1r,c_her2,c_erk)
trainingLabelihe<-list(control=list(igf1r=1:12,her2=1:12,erk=1:12),igf1r=13:18,her2=19:23,erk=24:29)
sub_dir=paste(basedir,"igf1r_her2_erk",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = trainihe,test=c_test,trainingLabel1 = trainingLabelihe,g=c(100,15,100),out_dir_base = sub_dir)

########################trying four at a time#####################
#1. AKT, BAD, HER2 & IGF1R
trainabhi<-cbind(c_gfp,c_akt,c_bad,c_her2,c_igf1r)
trainingLabelabhi<-list(control=list(akt=1:12,bad=1:12,her2=1:12,igf1r=1:12),akt=13:18, bad=19:24,her2=25:29,igf1r=30:35)
sub_dir=paste(basedir,"akt_bad_her2_igf1r",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = trainabhi,test=c_test,trainingLabel1 = trainingLabelabhi,g=c(150,150,15,100),out_dir_base = sub_dir)

#2. AKT, BAD, HER2 & ERH
trainabhe<-cbind(c_gfp,c_akt,c_bad,c_her2,c_erk)
trainingLabelabhe<-list(control=list(akt=1:12,bad=1:12,her2=1:12,erk=1:12),akt=13:18, bad=19:24,her2=25:29,erk=30:35)
sub_dir=paste(basedir,"akt_bad_her2_erk",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = trainabhe,test=c_test,trainingLabel1 = trainingLabelabhe,g=c(150,150,15,100),out_dir_base = sub_dir)

#3. AKT, BAD, IGF1R & ERK
trainabie<-cbind(c_gfp,c_akt,c_bad,c_igf1r,c_erk)
trainingLabelabie<-list(control=list(akt=1:12,bad=1:12,igf1r=1:12,erk=1:12),akt=13:18, bad=19:24,igf1r=25:30,erk=31:36)
sub_dir=paste(basedir,"akt_bad_igf1r_erk",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = trainabie,test=c_test,trainingLabel1 = trainingLabelabie,g=c(150,150,100,100),out_dir_base = sub_dir)

#4. AKT, IGF1R, HER2 & ERK
trainaihe<-cbind(c_gfp,c_akt,c_igf1r,c_her2,c_erk)
trainingLabelaihe<-list(control=list(akt=1:12,igf1r=1:12,her2=1:12,erk=1:12),akt=13:18, 1gf1r=19:24,her2=25:29,erk=30:35)
sub_dir=paste(basedir,"akt_igf1r_her2_erk",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = trainaihe,test=c_test,trainingLabel1 = trainingLabelaihe,g=c(150,100,15,100),out_dir_base = sub_dir)

#5. BAD, IGF1R, HER2 & ERK
trainbihe<-cbind(c_gfp,c_bad,c_igf1r,c_her2,c_erk)
trainingLabelbihe<-list(control=list(bad=1:12,igf1r=1:12,her2=1:12,erk=1:12),bad=13:18, 1gf1r=19:24,her2=25:29,erk=30:35)
sub_dir=paste(basedir,"bad_igf1r_her2_erk",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = trainbihe,test=c_test,trainingLabel1 = trainingLabelbihe,g=c(150,100,15,100),out_dir_base = sub_dir)


#########including all 5 pathways######
trainhall5<-cbind(c_gfp,c_akt,c_bad,c_her2,c_igf1r, c_erk)
trainingLabelall5<-list(control=list(akt=1:12,bad=1:12,her2=1:12,igf1r=1:12, erk=1:12),akt=13:18, bad=19:24,her2=25:29,igf1r=30:35, erk=36:41)
sub_dir=paste(basedir,"akt_bad_her2_igf1r_erk",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = trainhall5,test=c_test,trainingLabel1 = trainingLabelall5,g=c(150,150,15,100,100),out_dir_base = sub_dir)