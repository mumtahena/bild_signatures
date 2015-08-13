Now,combining predictions and proteomics data

data_icbp<-gatherFile("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/")

colnames(data_icbp)<-gsub(pattern = "/pathway_activity_testset.csv",replacement = "",x = colnames(data_icbp))
rownames(data_icbp)[1:7]<-c("184A1","184B5","21MT1","21MT2","21NT","21PT","600MPE")

prot<-read.table("~/Dropbox/bild_signatures/bild_signatures/Datasets/proteomics.txt",sep='\t',header=1,row.names=1)
pred_prot<-merge_drop(data_icbp,prot)
write.table(pred_prot,"~/Dropbox/bild_signatures/bild_signatures/Datasets/pred_130_prot_70.txt",sep='\t',col.names=NA,quote=F)
write.table(pred_prot,"~/Dropbox/bild_signatures/bild_signatures/Datasets/pred_380_prot_70.txt",sep='\t',col.names=NA,quote=F)



cor_mat=p_mat=matrix(0,ncol(data_icbp),70)
rownames(cor_mat)=rownames(p_mat)=colnames(pred_prot)[1:ncol(data_icbp)]
colnames(cor_mat)=colnames(p_mat)=colnames(pred_prot)[(ncol(data_icbp)+1):ncol(pred_prot)]

for(i in 1:ncol(data_icbp)){
  for(j in 1:70){
    temp=cor.test(pred_prot[,i],pred_prot[,j+ncol(data_icbp)],use="pairwise",method="spearman")
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
write.table(cor_mat,"Multi_single_protein_cor_mat_4_30.txt",col.names = NA,quote=F,sep='\t')
write.table(p_mat,"Multi_single_protein_p_mat_4_30.txt",col.names = NA,quote=F,sep='\t')
write.table(cor_p_mat,"Multi_single_protein_cor_p_mat_4_30.txt",col.names = NA,quote=F,sep='\t')

#Identifying the best single and multipathway predictions

single_pathway_best<-data_icbp[,c("akt_75_gene_list/adap_adap_single/akt","her2_15_gene_list/adap_adap_single/her2","igf1r_75_gene_list/adap_adap_single/igf1r","bad_200_gene_list/adap_adap_single/bad","erk_250_gene_list/adap_adap_single/erk","raf_100_gene_list/adap_adap_single/raf","egfr_25_gene_list/adap_adap_single/egfr","krasgv_300_gene_list/adap_adap_single/sigProtein","krasqh_300_gene_list/adap_adap_single/sigProtein")]#,"kraswt_300_gene_list/adap_adap_single/sigProtein")]

multi_pathway_best<-data_icbp[,c("akt_bad_her2_raf/adap_adap_multi/akt","akt_bad_her2_erk/adap_adap_multi/bad","her2_igf1r/adap_adap_multi/her2","igf1r_raf_her2_erk/adap_adap_multi/igf1r","erk_bad/adap_adap_multi/erk","akt_bad_igf1r_raf/adapB_multi/raf","egfr_25_gene_list/adap_adap_single/egfr","krasgv_300_gene_list/adap_adap_single/sigProtein","krasqh_300_gene_list/adap_adap_single/sigProtein")]#,"kraswt_300_gene_list/adap_adap_single/sigProtein")]
colnames(single_pathway_best)=c("AKT","HER2","IGF1R","BAD","ERK","RAF","EGFR","KRASGV","KRASQH")#,"KRASWT")           
#heatmap.2(as.matrix(single_pathway_best),col = bluered,density.info = "none",trace="none",margins = c(8,8),scale = "row")
colnames(multi_pathway_best)=c("AKT","BAD","HER2","IGF1R","ERK","RAF","EGFR","KRASGV","KRASQH")#,"KRASWT")           
#prot<-t(read.table("~/Dropbox/bild_signatures/Datasets/dream7_prot.txt",sep='\t',header=1,row.names=1,check.names = T))
#colnames(prot)<-gsub("-","",colnames(prot))

single_prot<-merge_drop(single_pathway_best,prot)
multi_prot<-merge_drop(multi_pathway_best,prot)

#---------------
par(mfrow=c(1,1),lwd = 5)

pdf("~/Dropbox/bild_signatures/protein_single_validation.pdf")
par(mfrow=c(1,1),lwd = 4)
tmp=cor.test(single_prot$AKT,single_prot$Akt,method="spearman",use="pairwise")
plot(single_prot$AKT,single_prot$Akt,ylab = "RPPA AKT Protein Score",xlab="Estimated AKT pathway activity",main=paste("Pearson's r=",round(tmp$estimate,digits = 2),"\np-value=",round(tmp$p.value,digits = 3),sep=" "))
abline(lm(single_prot$Akt~single_prot$AKT),col="red")

tmp=cor.test(single_prot$BAD,single_prot$Akt)
plot(single_prot$BAD,single_prot$Akt,ylab = "RPPA AKT Protein Score",xlab="Estimated BAD pathway activity",main=paste("Pearson's r=",round(tmp$estimate,digits = 2),"\np-value=",round(tmp$p.value,digits = 3),sep=" "))
abline(lm(single_prot$Akt~single_prot$BAD),col="red")

tmp=cor.test(single_prot$EGFR.x,single_prot$EGFRp1068)
plot(single_prot$EGFR.x,single_prot$EGFRp1068,ylab = "RPPA EGFRp1068 Protein Score",xlab="Estimated EGFR pathway activity",main=paste("Pearson's r=",round(tmp$estimate,digits = 2),"\np-value=",round(tmp$p.value,digits = 3),sep=" "))
abline(lm(single_prot$EGFRp1068~single_prot$EGFR.x),col="red")

tmp=cor.test(single_prot$ERK,single_prot$PKCalphap657)
plot(single_prot$ERK,single_prot$PKCalphap657,ylab = "RPPA PKCalphap657 Protein Score",xlab="Estimated ERK pathway activity",main=paste("Pearson's r=",round(tmp$estimate,digits = 2),"\np-value=",round(tmp$p.value,digits = 4),sep=" "))
abline(lm(single_prot$PKCalphap657~single_prot$ERK),col="red")

tmp=cor.test(single_prot$HER2.x,single_prot$HER2p1248)
plot(single_prot$HER2.x,single_prot$HER2p1248,ylab = "RPPA HER2p1248 Protein Score",xlab="Estimated HER2 pathway activity",main=paste("Pearson's r=",round(tmp$estimate,digits = 2),"\np-value=",round(tmp$p.value,digits = 7),sep=" "))
abline(lm(single_prot$HER2p1248~single_prot$HER2.x),col="red")

tmp=cor.test(single_prot$IGF1R,single_prot$IGFR1)
plot(single_prot$IGF1R,single_prot$IGFR1,ylab = "RPPA IGF1R Protein Score",xlab="Estimated IGF1R pathway activity",main=paste("Pearson's r=",round(tmp$estimate,digits = 2),"\np-value=",round(tmp$p.value,digits = 2),sep=" "))
abline(lm(single_prot$IGFR1~single_prot$IGF1R),col="red")

tmp=cor.test(single_prot$KRASGV,single_prot$EGFRp1068)
plot(single_prot$KRASGV,single_prot$EGFRp1068,ylab = "RPPA EGFRp1068 Protein Score",xlab="Estimated KRASG12V(mutant) pathway activity",main=paste("Pearson's r=",round(tmp$estimate,digits = 2),"\np-value=",round(tmp$p.value,digits = 2),sep=" "))
abline(lm(single_prot$EGFRp1068~single_prot$KRASGV),col="red")

tmp=cor.test(single_prot$KRASQH,single_prot$EGFRp1068)
plot(single_prot$KRASQH,single_prot$EGFRp1068,ylab = "RPPA EGFRp1068 Protein Score",xlab="Estimated KRASQ61H(mutant) pathway activity",main=paste("Pearson's r=",round(tmp$estimate,digits = 2),"\np-value=",round(tmp$p.value,digits = 2),sep=" "))
abline(lm(single_prot$EGFRp1068~single_prot$KRASQH),col="red")

tmp=cor.test(single_prot$MAPKp,single_prot$RAF)
plot(single_prot$RAF,single_prot$MAPKp,ylab = "RPPA MAPKp Protein Score",xlab="Estimated RAF pathway activity",main=paste("Pearson's r=",round(tmp$estimate,digits = 2),"\np-value=",round(tmp$p.value,digits = 2),sep=" "))
abline(lm(single_prot$MAPKp~single_prot$RAF),col="red")
dev.off()


pdf("~/Dropbox/protein_multi_validation.pdf")
par(lwd=5)
tmp=cor.test(multi_prot$AKT,multi_prot$Akt,method="spearman",use="pairwise")
plot(multi_prot$AKT,multi_prot$Akt,ylab = "RPPA AKT Protein Score",xlab="Estimated AKT pathway activity",main=paste("R=",round(tmp$estimate,digits = 2),"\np-value=",round(tmp$p.value,digits = 3),sep=" "))
abline(lm(multi_prot$Akt~multi_prot$AKT),col="red")

tmp=cor.test(multi_prot$BAD,multi_prot$Akt,method="spearman",use="pairwise")
plot(multi_prot$BAD,multi_prot$Akt,ylab = "RPPA AKT Protein Score",xlab="Estimated BAD pathway activity",main=paste("R=",round(tmp$estimate,digits = 2),"\np-value=",round(tmp$p.value,digits = 3),sep=" "))
abline(lm(multi_prot$Akt~multi_prot$BAD),col="red")

tmp=cor.test(multi_prot$EGFR.x,multi_prot$EGFRp1068,method="spearman",use="pairwise")
plot(multi_prot$EGFR.x,multi_prot$EGFRp1068,ylab = "RPPA EGFRp1068 Protein Score",xlab="Estimated EGFR pathway activity",main=paste("R=",round(tmp$estimate,digits = 2),"\np-value=",round(tmp$p.value,digits = 3),sep=" "))
abline(lm(multi_prot$EGFRp1068~multi_prot$EGFR.x),col="red")

tmp=cor.test(multi_prot$ERK,multi_prot$PKCalphap657,method="spearman",use="pairwise")
plot(multi_prot$ERK,multi_prot$PKCalphap657,ylab = "RPPA PKCalphap657 Protein Score",xlab="Estimated ERK pathway activity",main=paste("R=",round(tmp$estimate,digits = 2),"\np-value=",round(tmp$p.value,digits = 4),sep=" "))
abline(lm(multi_prot$PKCalphap657~multi_prot$ERK),col="red")

tmp=cor.test(multi_prot$HER2.x,multi_prot$HER2p1248,method="spearman",use="pairwise")
plot(multi_prot$HER2.x,multi_prot$HER2p1248,ylab = "RPPA HER2p1248 Protein Score",xlab="Estimated HER2 pathway activity",main=paste("R=",round(tmp$estimate,digits = 2),"\np-value=",round(tmp$p.value,digits = 7),sep=" "))
abline(lm(multi_prot$HER2p1248~multi_prot$HER2.x),col="red")

tmp=cor.test(multi_prot$IGF1R,multi_prot$PDK1p241,method="spearman",use="pairwise")
plot(multi_prot$IGF1R,multi_prot$PDK1p241,ylab = "RPPA PDK1p241 Protein Score",xlab="Estimated IGF1R pathway activity",main=paste("R=",round(tmp$estimate,digits = 2),"\np-value=",round(tmp$p.value,digits = 2),sep=" "))
abline(lm(multi_prot$IGFR1~multi_prot$IGF1R),col="red")

tmp=cor.test(multi_prot$KRASGV,multi_prot$EGFR.y,method="spearman",use="pairwise")
plot(multi_prot$KRASGV,multi_prot$EGFR.y,ylab = "RPPA EGFRp1068 Protein Score",xlab="Estimated KRASG12V(mutant) pathway activity",main=paste("R=",round(tmp$estimate,digits = 2),"\np-value=",round(tmp$p.value,digits = 2),sep=" "))
abline(lm(multi_prot$EGFR.y~multi_prot$KRASGV),col="red")

tmp=cor.test(multi_prot$KRASQH,multi_prot$EGFR.y,method="spearman",use="pairwise")
plot(multi_prot$KRASQH,multi_prot$EGFR.y,ylab = "RPPA EGFRp1068 Protein Score",xlab="Estimated KRASQ61H(mutant) pathway activity",main=paste("R=",round(tmp$estimate,digits = 2),"\np-value=",round(tmp$p.value,digits = 2),sep=" "))
abline(lm(multi_prot$EGFR.y~multi_prot$KRASQH),col="red")

tmp=cor.test(multi_prot$MAPKp,multi_prot$RAF,method="spearman",use="pairwise")
plot(multi_prot$RAF,multi_prot$MAPKp,ylab = "RPPA MAPKp Protein Score",xlab="Estimated RAF pathway activity",main=paste("R=",round(tmp$estimate,digits = 2),"\np-value=",round(tmp$p.value,digits = 2),sep=" "))
abline(lm(multi_prot$MAPKp~multi_prot$RAF),col="red")
dev.off()


data_tcga_br<-gatherFile("~/Desktop/PANCAN24_BRCA_1119/")
colnames(data_tcga_br)<-gsub(pattern = "/pathway_activity_testset.csv",replacement = "",x = colnames(data_tcga_br))
tcga_single_pathway_best<-data_tcga_br[,c("akt_75_gene_list/adap_adap_single/akt","her2_15_gene_list/adap_adap_single/her2","igf1r_75_gene_list/adap_adap_single/igf1r","bad_200_gene_list/adap_adap_single/bad","erk_250_gene_list/adap_adap_single/erk","raf_100_gene_list/adap_adap_single/raf","egfr_25_gene_list/adap_adap_single/egfr","krasgv_300_gene_list/adap_adap_single/sigProtein","krasqh_300_gene_list/adap_adap_single/sigProtein")]#,"kraswt_300_gene_list/adap_adap_single/sigProtein")]

multi_pathway_best<-data_icbp[,c("akt_bad_her2_raf/adap_adap_multi/akt","akt_bad_her2_erk/adap_adap_multi/bad","her2_igf1r/adap_adap_multi/her2","igf1r_raf_her2_erk/adap_adap_multi/igf1r","erk_bad/adap_adap_multi/erk","akt_bad_igf1r_raf/adapB_multi/raf","egfr_25_gene_list/adap_adap_single/egfr","krasgv_300_gene_list/adap_adap_single/sigProtein","krasqh_300_gene_list/adap_adap_single/sigProtein")]#,"kraswt_300_gene_list/adap_adap_single/sigProtein")]
colnames(single_pathway_best)=c("AKT","HER2","IGF1R","BAD","ERK","RAF","EGFR","KRASGV","KRASQH")#,"KRASWT")           
#heatmap.2(as.matrix(single_pathway_best),col = bluered,density.info = "none",trace="none",margins = c(8,8),scale = "row")
colnames(multi_pathway_best)=c("AKT","BAD","HER2","IGF1R","ERK","RAF","EGFR","KRASGV","KRASQH")#,"KRASWT")           
#prot<-t(read.table("~/Dropbox/bild_signatures/Datasets/dream7_prot.txt",sep='\t',header=1,row.names=1,check.names = T))
#colnames(prot)<-gsub("-","",colnames(prot))

par(lwd=2)
heatmap.2(as.matrix(cor(multi_pathway_best,method="spearman",use="pairwise")),col = bluered,density.info = "none",trace="none",margins = c(8,8),scale = "row")

heatmap.2(as.matrix(cor(multi_pathway_best_tcga_wo_ERK,method="spearman",use="pairwise")),col = bluered,density.info = "none",trace="none",margins = c(8,8),scale = "row")





