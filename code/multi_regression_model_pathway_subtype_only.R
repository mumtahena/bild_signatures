setwd("~/Desktop/drug reponse analysis/")
#dat = read.table("icbp_preds_drug_combat_adap.txt",header=T,row.names=1,sep='\t')
multi<-read.table("~/Desktop/multipathway_preds", sep='\t',row.names=1, header=T)
single<-read.csv("~/Desktop/tmp/single_pathway_results.csv", row.names=1, header=T)
dat<-cbind(multi,single)
rownames(dat)[1:7]<-c("184A1","184B5","21MT1","21MT2","21NT","21PT","600MPE")
drugs<-read.delim("~/Dropbox/bild_signatures/Datasets/ICBP_drugs.txt", header=1, sep='\t',row.names=1)
merge_drop<-function(x,y,by=0)
{
  new_m<-merge(x,y,by=by)
  rownames(new_m)<-new_m$Row.names
  return(new_m[,2:length(colnames(new_m))])
}


dat<-merge_drop(drugs,dat,by=0)
subtype = dat[,1]

colnames(dat)
snp <- t(read.table("snp_all_v2.txt", row.names=1, header=T)) #71 250
indel <- read.table("indel_all.txt", row.names=1,header=T) #36 38

for (i in 1:nrow(indel)){
  for (j in 1:ncol(indel)){
    iindex = which(rownames(snp)==rownames(indel)[i])
    if (!(colnames(indel)[j] %in% colnames(snp))){snp=cbind(snp,0);colnames(snp)[ncol(snp)]=colnames(indel)[j]}
    jindex = which(colnames(snp)==colnames(indel)[j])
    snp[iindex,jindex]=min(snp[iindex,jindex]+indel[i,j],1)
  }
}



proteomics = read.table("proteomics.txt", row.names=1,header=T)


cell = rownames(dat)
subtype = dat[,1]
subtype_ERBB2 = dat[,2]


snps=snp[match(cell,rownames(snp)),]
prot=proteomics[match(cell,rownames(proteomics)),]


#pathway = dat[,93:99]
#colnames(pathway)=c('her2','raf','bad','akt','akt5','igf1r','erk')
pathway = dat[,101:148]
#colnames(pathway)=c('akt_multi','bad_multi','her2_multi','erk_multi','igf1r_multi','raf_multi')

#### Pathway by subtype ####
# pdf(file='activity_subtype.pdf',width=15,height=7.5)
# par(mfrow=c(2,1))
# for (i in 1:ncol(pathway)){
# 	plot(pathway[,i]~subtype)
# 	title(colnames(pathway)[i])
# }
# dev.off()
# 
# pdf(file='activity_subtype_ERBB2.pdf',width=15,height=7.5)
# par(mfrow=c(2,1))
# for (i in 1:ncol(pathway)){
# 	plot(pathway[,i]~subtype_ERBB2)
# 	title(colnames(pathway)[i])
# }
# dev.off()


#### Corellation analysis ####
drugs=dat[,11:100]
cors = matrix(0,ncol(drugs),ncol(pathway))
rownames(cors)=colnames(drugs)
colnames(cors)=colnames(pathway)

for (i in 1:ncol(pathway)){
  for (j in 1:ncol(drugs)){
    cors[j,i]=cor(pathway[,i],drugs[,j],use="pairwise.complete.obs")
  }
}

cors
write.table(cors,"~/Desktop/all_cors.txt", sep='\t',col.names=NA,quote=F)

### Regression analysis ###
multipathway = dat[,137:140]
singlepathway= dat[,145:148]
preds<-dat[,101:148]
##colnames(pathway)=c('akt_multi','bad_multi','her2_multi','erk_multi','igf1r_multi','raf_multi','akt_single','bad_single','her2_single','erk_single','igf1r_single','raf_single')
#pathway = pathway[,1:6]
#pathway = pathway[,7:12]

library(MASS)
drugs=dat[,11:100]
#shortlist= c('Everolimus', "Erlotinib")
shortlist= c('Sigma.AKT1.2.inhibitor','BEZ235','BIBW2992','Everolimus','GSK2119563','GSK2126458','GSK2141795','GSK1059615','GSK650394','Lapatinib')
shortlist= c('Sigma.AKT1.2.inhibitor','GSK2141795','Triciribine','Everolimus','Temsirolimus','Lapatinib','BEZ235','BIBW2992')
shortlist=c("Sigma.AKT1.2.inhibitor",'GSK2141795','Triciribine','Everolimus','Temsirolimus','Lapatinib')
#shortlist = colnames(drugs)
drugs_short=drugs[,match(shortlist,colnames(drugs))]
#rownames(drugs_short)=cell

#### results ####
library(MASS)
############for multipathway predictions########
multiresults=multimodels=list()
for (i in 1:ncol(drugs_short)){
  fit1 <- lm(drugs_short[,i]~.,data=multipathway)
  step1 <- stepAIC(fit1, direction="both")
  
  fit2 <- lm(drugs_short[,i]~.+subtype,data=multipathway)
  step2 <- stepAIC(fit2, direction="both")
  
  fit3 <- lm(drugs_short[,i]~subtype)
  
  multiresults[[shortlist[i]]]=list(R2=summary(step1)$r.squared,model=round(summary(step1)$coeff,4),R2_sub=summary(step2)$r.squared,model_sub=round(summary(step2)$coeff,4),R2_sub_only=summary(fit3)$r.squared,model_sub_only=round(summary(fit3)$coeff,4))
  multimodels[[shortlist[i]]]=list(mod=step1,mod_sub=step2,mod_sub_only=fit3)
}
setwd("~/Dropbox/bild_signatures/bild_signatures/Results/")
write(file="multipathway_results.xls","Results for MultiPathway predictions/Drug Models")
for (i in 1:length(multiresults)){
  write('',file="multipathway_results.xls", append=T)
  write(names(multiresults)[i], file="multipathway_results.xls", append=T)
  write(paste("Model R2 (no subtype):",round(multiresults[[i]]$R2,4), sep=' '),file="multipathway_results.xls", append=T)
  cat('\t',file="multipathway_results.xls", append=T)
  write.table(multiresults[[i]]$model,file="multipathway_results.xls", quote=F, col.names=T, row.names=T, append=T,sep='\t')
  write('',file="multipathway_results.xls", append=T)
  write(paste("Model R2 (+subtype):",round(multiresults[[i]]$R2_sub,4), sep=' '),file="multipathway_results.xls", append=T)
  cat('\t',file="multipathway_results.xls", append=T)
  write.table(multiresults[[i]]$model_sub,file="multipathway_results.xls", quote=F, col.names=T, row.names=T, append=T,sep='\t')
  write('',file="multipathway_results.xls", append=T)
  write(paste("Model subtype only:",round(multiresults[[i]]$R2_sub_only,4), sep=' '),file="multipathway_results.xls", append=T)
  cat('\t',file="multipathway_results.xls", append=T)
  write.table(multiresults[[i]]$model_sub_only,file="multipathway_results.xls", quote=F, col.names=T, row.names=T, append=T,sep='\t')
  
}
###########for single pathway predictions#############
singleresults=singlemodels=list()
for (i in 1:ncol(drugs_short)){
  fit1 <- lm(drugs_short[,i]~.,data=singlepathway)
  step1 <- stepAIC(fit1, direction="both")
  
  fit2 <- lm(drugs_short[,i]~.+subtype,data=singlepathway)
  step2 <- stepAIC(fit2, direction="both")
  
  fit3 <- lm(drugs_short[,i]~subtype)
  
  singleresults[[shortlist[i]]]=list(R2=summary(step1)$r.squared,model=round(summary(step1)$coeff,4),R2_sub=summary(step2)$r.squared,model_sub=round(summary(step2)$coeff,4),R2_sub_only=summary(fit3)$r.squared,model_sub_only=round(summary(fit3)$coeff,4))
  singlemodels[[shortlist[i]]]=list(mod=step1,mod_sub=step2,mod_sub_only=fit3)
}

write(file="Single pathway_results.xls","Results for Single Pathway predictions/Drug Models")
for (i in 1:length(singleresults)){
  write('',file="singlepathway_results.xls", append=T)
  write(names(singleresults)[i], file="singlepathway_results.xls", append=T)
  write(paste("Model R2 (no subtype):",round(singleresults[[i]]$R2,4), sep=' '),file="singlepathway_results.xls", append=T)
  cat('\t',file="singlepathway_results.xls", append=T)
  write.table(singleresults[[i]]$model,file="singlepathway_results.xls", quote=F, col.names=T, row.names=T, append=T,sep='\t')
  write('',file="singlepathway_results.xls", append=T)
  write(paste("Model R2 (+subtype):",round(singleresults[[i]]$R2_sub,4), sep=' '),file="singlepathway_results.xls", append=T)
  cat('\t',file="singlepathway_results.xls", append=T)
  write.table(singleresults[[i]]$model_sub,file="singlepathway_results.xls", quote=F, col.names=T, row.names=T, append=T,sep='\t')
  
}
