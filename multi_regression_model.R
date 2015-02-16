setwd("~/Desktop/drug reponse analysis/")
#dat = read.table("icbp_preds_drug_combat_adap.txt",header=T,row.names=1,sep='\t')
multi<-read.csv("~/Desktop/tmp/multi_pathway_results.csv", row.names=1, header=T)
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
pathway = dat[,101:104]
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


### Regression analysis ###
multipathway = dat[,101:104]
singlepathway= dat[,105:108]
##colnames(pathway)=c('akt_multi','bad_multi','her2_multi','erk_multi','igf1r_multi','raf_multi','akt_single','bad_single','her2_single','erk_single','igf1r_single','raf_single')
#pathway = pathway[,1:6]
#pathway = pathway[,7:12]

library(MASS)
drugs=dat[,11:100]
#shortlist= c('Everolimus', "Erlotinib")
#shortlist= c('Sigma.AKT1.2.inhibitor','BEZ235','BIBW2992','Everolimus','GSK2119563','GSK2126458','GSK2141795','GSK1059615','GSK650394','Lapatinib')
shortlist= c('Sigma.AKT1.2.inhibitor','GSK2141795','Triciribine')

#shortlist = colnames(drugs)
drugs_short=drugs[,match(shortlist,colnames(drugs))]
rownames(drugs_short)=cell

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



#### Now for the snps!
r2_cut = 0.15
r2_add_cut = 0.15

for (j in 1:ncol(drugs_short)){
  mod=models[[j]]$mod
  keep = cell %in% rownames(mod$model)
  snp_short = snps[keep,]
  cat(paste("Using",dim(snp_short)[1], 'samples from a total of', dim(snps)[1],'\n'))
  
  snp_keep = apply(snp_short,2,sum,na.rm=T)>3
  snp_shorter = snp_short[,snp_keep]
  cat(paste("Using",dim(snp_shorter)[2], 'SNPs from a total of', dim(snps)[2],'\n'))
  print(colnames(snp_shorter))
  
  r2 = r2_alone = r2_added = NULL
  for (i in 1:ncol(snp_shorter)){
    r2=c(r2, summary(mod)$r.squared)
    r2_alone=c(r2_alone, summary(lm(mod$model[,1]~snp_shorter[,i]))$r.squared)
    r2_added=c(r2_added, summary(lm(mod$model[,1]~.+snp_shorter[,i],data=data.frame(1,mod$model[,-1])))$r.squared)
  }
  names(r2_alone)= names(r2_alone)= names(r2_added)= colnames(snp_shorter)
  r2_alone = sort(r2_alone,T)
  r2_diff = sort(r2_added-r2,T)
  
  results[[shortlist[j]]]$genes=paste(names(which(r2_alone>r2_cut)),' (',round(r2_alone,2)[which((r2_alone)>r2_cut)],'),',sep='')
  results[[shortlist[j]]]$genes_add=paste(names(which(r2_diff>r2_add_cut)),' (',round(r2_diff,2)[which((r2_diff)>r2_add_cut)],'),',sep='')
  models[[shortlist[j]]]$genes_best = names(r2_diff)[1]
  
  mod=models[[j]]$mod_sub
  
  r2 = r2_added = NULL
  for (i in 1:ncol(snp_shorter)){
    r2=c(r2, summary(mod)$r.squared)
    r2_added=c(r2_added, summary(lm(mod$model[,1]~.+snp_shorter[,i],data=data.frame(1,mod$model[,-1])))$r.squared)
  }
  names(r2)= names(r2_added)= colnames(snp_shorter)
  r2_diff = sort(r2_added-r2,T)
  
  results[[shortlist[j]]]$genes_add_sub=paste(names(which(r2_diff>r2_add_cut)),' (',round(r2_diff,2)[which((r2_diff)>r2_add_cut)],'),',sep='')
  models[[shortlist[j]]]$genes_best_sub = names(r2_diff)[1]
}


#### RPPA  #####
r2_cut = 0.15
r2_add_cut = 0.15

for (j in 1:ncol(drugs_short)){
  mod=models[[j]]$mod
  keep = cell %in% rownames(mod$model)
  prot_short = prot[keep,]
  cat(paste("Using",dim(prot_short)[1], 'samples from a total of', dim(proteomics)[1],'\n'))
  
  r2 = r2_alone = r2_added = NULL
  for (i in 1:ncol(prot_short)){
    r2=c(r2, summary(mod)$r.squared)
    r2_alone=c(r2_alone, summary(lm(mod$model[,1]~prot_short[,i]))$r.squared)
    r2_added=c(r2_added, summary(lm(mod$model[,1]~.+prot_short[,i],data=data.frame(1,mod$model[,-1])))$r.squared)
  }
  names(r2_alone)= names(r2_alone)= names(r2_added)= colnames(prot_short)
  r2_alone = sort(r2_alone,T)
  r2_diff = sort(r2_added-r2,T)
  
  results[[shortlist[j]]]$proteins=paste(names(which(r2_alone>r2_cut)),' (',round(r2_alone,2)[which((r2_alone)>r2_cut)],'),',sep='')
  results[[shortlist[j]]]$proteins_add=paste(names(which(r2_diff>r2_add_cut)),' (',round(r2_diff,2)[which((r2_diff)>r2_add_cut)],'),',sep='')
  models[[shortlist[j]]]$proteins_best = names(r2_diff)[1]
  
  mod=models[[j]]$mod_sub
  r2 = r2_added = NULL
  for (i in 1:ncol(prot_short)){
    r2=c(r2, summary(mod)$r.squared)
    r2_added=c(r2_added, summary(lm(mod$model[,1]~.+prot_short[,i],data=data.frame(1,mod$model[,-1])))$r.squared)
  }
  names(r2)= names(r2_added)= colnames(prot_short)
  r2_diff = sort(r2_added-r2,T)
  
  results[[shortlist[j]]]$proteins_add_sub=paste(names(which(r2_diff>r2_add_cut)),' (',round(r2_diff,2)[which((r2_diff)>r2_add_cut)],'),',sep='')
  models[[shortlist[j]]]$proteins_best_sub = names(r2_diff)[1]
}




### Write results:

write(file="results.xls","Results for Pathway/Drug Models")
for (i in 1:length(results)){
  write('',file="results.xls", append=T)
  write(names(results)[i], file="results.xls", append=T)
  write(paste("Model R2 (no subtype):",round(results[[i]]$R2,4), sep=' '),file="results.xls", append=T)
  cat('\t',file="results.xls", append=T)
  write.table(results[[i]]$model,file="results.xls", quote=F, col.names=T, row.names=T, append=T,sep='\t')
  write(paste(c("Important genes:",results[[i]]$genes), sep=' ',collapse=' '),file="results.xls", append=T)
  write(paste(c("Genes that add to the pathway model:",results[[i]]$genes_add), sep=' ',collapse=' '),file="results.xls", append=T)
  write(paste(c("Important proteins:",results[[i]]$proteins), sep=' ',collapse=' '),file="results.xls", append=T)
  write(paste(c("Proteins that add to the pathway model:",results[[i]]$proteins_add), sep=' ',collapse=' '),file="results.xls", append=T)
  
  write('',file="results.xls", append=T)
  write(paste("Model R2 (+subtype):",round(results[[i]]$R2_sub,4), sep=' '),file="results.xls", append=T)
  cat('\t',file="results.xls", append=T)
  write.table(results[[i]]$model_sub,file="results.xls", quote=F, col.names=T, row.names=T, append=T,sep='\t')
  write(paste(c("Important genes:",results[[i]]$genes), sep=' ',collapse=' '),file="results.xls", append=T)
  write(paste(c("Genes that add to the pathway model:",results[[i]]$genes_add_sub), sep=' ',collapse=' '),file="results.xls", append=T)
  write(paste(c("Important proteins:",results[[i]]$proteins), sep=' ',collapse=' '),file="results.xls", append=T)
  write(paste(c("Proteins that add to the pathway model:",results[[i]]$proteins_add_sub), sep=' ',collapse=' '),file="results.xls", append=T)
  write('',file="results.xls", append=T)
}



#### Plots

par(mfrow = c(1,1))
plot(models$'Sigma.AKT1.2.inhibitor'$mod$fit,models$'Sigma.AKT1.2.inhibitor'$mod$model[,1],main=paste("Sigma.AKT1.2.inhibitor Predicted vs Actual Efficacy (R2=",round(summary(models$'Sigma.AKT1.2.inhibitor'$mod)$r.squared,2),")",sep=''),xlab="Precticted Efficacy",ylab="Actual Efficacy")
abline(0,1,col=2)
plot(models$'Sigma.AKT1.2.inhibitor'$mod_sub$fit,models$'Sigma.AKT1.2.inhibitor'$mod_sub$model[,1],main=paste("Sigma.AKT1.2.inhibitor Predicted vs Actual Efficacy (R2=",round(summary(models$'Sigma.AKT1.2.inhibitor'$mod_sub)$r.squared,2),")",sep=''),xlab="Precticted Efficacy",ylab="Actual Efficacy")
abline(0,1,col=2)

plot(models$'GSK2141795'$mod$fit,models$'GSK2141795'$mod$model[,1],main=paste("GSK2141795 Predicted vs Actual Efficacy (R2=",round(summary(models$'GSK2141795'$mod)$r.squared,2),")",sep=''),xlab="Precticted Efficacy",ylab="Actual Efficacy")
abline(0,1,col=2)
plot(models$'GSK2141795'$mod_sub$fit,models$'GSK2141795'$mod_sub$model[,1],main=paste("GSK2141795 Predicted vs Actual Efficacy (R2=",round(summary(models$'GSK2141795'$mod_sub)$r.squared,2),")",sep=''),xlab="Precticted Efficacy",ylab="Actual Efficacy")
abline(0,1,col=2)


plot(models$'Triciribine'$mod$fit,models$'Triciribine'$mod$model[,1],main=paste("Triciribine Predicted vs Actual Efficacy (R2=",round(summary(models$'Triciribine'$mod)$r.squared,2),")",sep=''),xlab="Precticted Efficacy",ylab="Actual Efficacy")
abline(0,1,col=2)
plot(models$'Triciribine'$mod_sub$fit,models$'Triciribine'$mod_sub$model[,1],main=paste("Triciribine Predicted vs Actual Efficacy (R2=",round(summary(models$'Triciribine'$mod_sub)$r.squared,2),")",sep=''),xlab="Precticted Efficacy",ylab="Actual Efficacy")
abline(0,1,col=2)

###################
plot(models$'Erlotinib'$mod$fit,models$'Erlotinib'$mod$model[,1],main=paste("Erlotinib Predicted vs Actual Efficacy (R2=",round(summary(models$'Erlotinib'$mod)$r.squared,2),")",sep=''),xlab="Precticted Efficacy",ylab="Actual Efficacy")
abline(0,1,col=2)
plot(models$'Erlotinib'$mod_sub$fit,models$'Erlotinib'$mod$model[,1],main=paste("Erlotinib (+subtype) Predicted vs Actual Efficacy (R2=",round(summary(models$'Erlotinib'$mod_sub)$r.squared,2),")",sep=''),xlab="Precticted Efficacy",ylab="Actual Efficacy")
abline(0,1,col=2)

plot(models$'Sorafenib'$mod$fit,models$'Sorafenib'$mod$model[,1],main=paste("Sorafenib Predicted vs Actual Efficacy (R2=",round(summary(models$'Sorafenib'$mod)$r.squared,2),")",sep=''),xlab="Precticted Efficacy",ylab="Actual Efficacy")
abline(0,1,col=2)
plot(models$'Sorafenib'$mod_sub$fit,models$'Sorafenib'$mod$model[,1],main=paste("Sorafenib (+subtype) Predicted vs Actual Efficacy (R2=",round(summary(models$'Sorafenib'$mod_sub)$r.squared,2),")",sep=''),xlab="Precticted Efficacy",ylab="Actual Efficacy")
abline(0,1,col=2)

#### Everolimus
mod=models$'Everolimus'$mod
keep = cell %in% rownames(mod$model)
prot_short = prot[keep,]

sub_keep=subtype[keep]
mod_sub=lm(mod$model[,1]~.+sub_keep, data=data.frame(1,mod$model[,-1]))


### add p21
mod_p21 = lm(mod$model[,1]~.+sub_keep+prot_short$p21, data=data.frame(1,mod$model[,-1]))

par(mfrow = c(1,3))
plot(models$'Everolimus'$mod$fit,models$'Everolimus'$mod$model[,1],main=paste("Everolimus Predicted vs Actual Efficacy (R^2=",round(summary(models$'Everolimus'$mod)$r.squared,2),")",sep=''),xlab="Precticted Efficacy",ylab="Actual Efficacy")
abline(0,1,col=2)
plot(mod_sub$fit,mod_sub$model[,1],main=paste("Everolimus (+subtype) Predicted vs Actual Efficacy (R2=",round(summary(mod_sub)$r.squared,2),")",sep=''),xlab="Precticted Efficacy",ylab="Actual Efficacy")
abline(0,1,col=2)
plot(mod_p21$fit,mod_p21$model[,1],main=paste("Everolimus (+subtype +p21) Predicted vs Actual Efficacy (R2=",round(summary(mod_p21)$r.squared,2),")",sep=''),xlab="Precticted Efficacy",ylab="Actual Efficacy")
abline(0,1,col=2)

cbind(mod$model[,1],mod$'fit')
cbind(mod_sub$model[,1],mod_sub$'fit')
cbind(mod_p21$model[,1],mod_p21$'fit')



#### more results
#
mod_GSK=models$'GSK2119563'$mod
mod_GSK_sub=models$'GSK2119563'$mod_sub

keep = cell %in% rownames(mod_GSK$model)
prot_short = prot[keep,]
mod_CCND1 = lm(mod_GSK_sub$model[,1]~.+prot_short$CCND1, data=data.frame(1,mod_GSK_sub$model[,-1]))

cbind(mod_GSK$model[,1],mod_GSK$'fit')
cbind(mod_GSK_sub$model[,1],mod_GSK_sub$'fit')
cbind(mod_CCND1$model[,1],mod_CCND1$'fit')

#
mod_AKT=models$'Sigma.AKT1.2'$mod
mod_AKT_sub=models$'Sigma.AKT1.2'$mod_sub

cbind(mod_AKT$model[,1],mod_AKT$'fit')
cbind(mod_AKT_sub$model[,1],mod_AKT_sub$'fit')



######  Pathway, snp, proteomics ....

snp_keep = apply(snps,2,sum,na.rm=T)>3
snp_short = snps[,snp_keep]

R2_pathway = numeric(ncol(drugs)) #matrix(0,ncol(drugs),ncol(pathway))
R2_subtype = numeric(ncol(drugs))
R2_protein = matrix(0,ncol(drugs),ncol(prot))
R2_snp = matrix(0,ncol(drugs),ncol(snp_short))

for (j in 1:ncol(drugs)){
  R2_subtype[j] = summary(lm(drugs[,j]~subtype))$r.squared
  
  fit_path = lm(drugs[,j]~.,data=pathway)
  step_path <- stepAIC(fit_path, direction="both")
  R2_pathway[j] = summary(step_path)$r.squared
  
  for (i in 1:ncol(pathway)){
    R2_snp[j,i] = summary(lm(drugs[,j]~snp_short[,i]))$r.squared
  }
  
  for (i in 1:ncol(prot)){
    R2_protein[j,i] = summary(lm(drugs[,j]~prot[,i]))$r.squared
  }
  
}

snp_max = apply(R2_snp,1,max)
prot_max = apply(R2_protein,1,max)

my_line <- function(x,y,...){
  points(x,y,...)
  abline(a = 0,b = 1,col=2)
}


pairs(cbind(R2_subtype,snp_max,R2_pathway,prot_max), lower.panel = my_line,upper.panel = my_line)

## Beat 
ind = which(colnames(drugs)=="Everolimus")
#25
bars = c(R2_pathway[25],R2_subtype[25],prot_max[25],snp_max[25])
names(bars) = c('Pathway',"Subtype",'Protein','Variant')

pdf('plots/everolimus1.pdf',7,7)
barplot(bars,ylab=expression(paste('Model ',R^2)), main=colnames(drugs)[ind],col=2:5)
dev.off()

pdf('plots/everolimus2.pdf',7,7)
barplot(bars[4:1],xlab=expression(paste('Model ',R^2)), main=colnames(drugs)[ind],horiz=T,col=5:2)
dev.off()

### combined 
full_mods = c(summary(mod_p21)$r.squared,summary(mod_CCND1)$r.squared,summary(mod_AKT_sub)$r.squared)
noprot_mods = c(summary(mod_sub)$r.squared,summary(mod_GSK_sub)$r.squared,summary(mod_AKT_sub)$r.squared)
path_mods = c(summary(mod)$r.squared,summary(mod_GSK)$r.squared,summary(mod_AKT)$r.squared)
names(noprot_mods)=names(full_mods)=names(path_mods)=c("Everolimus",'GSK2119563','Sigma.AKT1.2')

pdf('plots/several1.pdf',7,7)
barplot(full_mods,ylab=expression(paste('Model ',R^2)),col=4,ylim=c(0,1))
barplot(noprot_mods,ylab=expression(paste('Model ',R^2)),col=3,add=T)
barplot(path_mods,ylab=expression(paste('Model ',R^2)),add=T,col=2)
legend(2.5,1,legend=c('Pathway','Subtype','Protein'),col=2:4,lwd=5)
dev.off()

pdf('plots/several2.pdf',7,7)
barplot(full_mods,ylab=expression(paste('Model ',R^2)),col=4,horiz=T,xlim=c(0,1))
barplot(noprot_mods,ylab=expression(paste('Model ',R^2)),col=3,add=T,horiz=T)
barplot(path_mods,ylab=expression(paste('Model ',R^2)),add=T,col=2,horiz=T)
legend(.65,2.3,legend=c('Pathway','Subtype','Protein'),col=2:4,lwd=5)
dev.off()


### combined 

r2_path_sub_prot_snp = r2_path_sub_prot = r2_path_sub = r2_path = numeric()

for (j in 1:ncol(drugs)){
  mod = models[[j]]$mod
  snp_name = models[[j]]$genes_best
  prot_name = models[[j]]$proteins_best
  path_names = colnames(mod$model)[-1]
  
  drug_dat = drugs[,j] 
  snp_best = snps[,colnames(snps)==snp_name]
  prot_best = prot[,colnames(prot)==prot_name]
  path_best = as.matrix(cbind(1,pathway[,colnames(pathway) %in% path_names]))
  
  r2_path_sub_prot_snp = c(r2_path_sub_prot_snp, summary(lm(drug_dat~path_best+subtype+snp_best+prot_best))$r.squared)
  r2_path_sub_prot = c(r2_path_sub_prot, summary(lm(drug_dat~path_best+subtype+prot_best))$r.squared)
  r2_path_sub = c(r2_path_sub, summary(lm(drug_dat~path_best+subtype))$r.squared)
  r2_path = c(r2_path, summary(lm(drug_dat~path_best))$r.squared)
}

ord = order(apply(cbind(r2_path_sub_prot_snp,r2_path,r2_path_sub_prot,r2_path_sub),1,max),decreasing=T)
names(r2_path_sub_prot_snp) = colnames(drugs)


pdf('plots/all.pdf',14,7)
par(las=2) # make label text perpendicular to axis
barplot(r2_path_sub_prot_snp[ord],col=5,cex.names=.65,ylab=expression(paste('Model ',R^2)), ylim=c(0,1))
barplot(r2_path_sub_prot[ord],col=4,add=T)
barplot(r2_path_sub[ord],col=3,add=T)
barplot(r2_path[ord],col=2,add=T)
legend(50,1,legend=c('Pathway','Subtype','Protein',"Variant"),col=2:5,lwd=5)
dev.off()

pdf('plots/all2.pdf',7,14)
par(las=2) # make label text perpendicular to axis
par(mar=c(5,8,4,2)) # increase y-axis margin.
barplot(sumr2[order(sumr2)],horiz=T,col=5,cex.names=.5,xlab=expression(paste('Model ',R^2)))
barplot((R2_subtype/2+R2_pathway+prot_max/2)[order(sumr2)],col=4,add=T,horiz=T)
barplot((R2_subtype/2+R2_pathway)[order(sumr2)],col=3,add=T,horiz=T)
barplot((R2_pathway)[order(sumr2)],col=2,add=T,horiz=T)
legend(.75,50,legend=c('Pathway','Subtype','Protein',"Variant"),col=2:5,lwd=5)
dev.off()
