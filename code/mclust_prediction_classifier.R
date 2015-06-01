drugs_srt<-drugs[rownames(drugs)%in%com$Cell.lines,c(15,49,33,77,27,26,32,40)]
library(mclust)
library(reshape2)
library(ggplot2)
akt_preds<-subset(pred_drug,select=grep("/akt", colnames(pred_drug) , ignore.case=FALSE, fixed=T))
bad_preds<-subset(pred_drug,select=grep("/bad", colnames(pred_drug) , ignore.case=FALSE, fixed=T))
labels=gsub("/adapB_multi/akt","AB",colnames(akt_preds))
labels=gsub("/adap_adap_multi/akt","AA",labels)
head(melt(akt_preds,labels="Labels"))
pdf("test_temp_akt.pdf")
boxplot(pred_drug$Sigma.AKT1.2.inhibitor~classes$classification,names = labels[1],col=c("red","green"),horizontal = T)

apply(akt_preds,2,Mclust(G=1:2))
for(i in 1:3){#ncol(akt_preds)){
  classes=Mclust(akt_preds[,i],G=1:2)
  #print(classes$classification)
  if(classes$G==2)
    {
    
    barplot(pred_drug$Sigma.AKT1.2.inhibitor~classes$classification,names=paste(labels[i],c("Low","High"),sep=":"),notch=T,col=c("red","green"),horizontal = T,type='l',add=T)
    #mtext(side = 3,text = paste("p-val:",round(t.test(pred_drug$Sigma.AKT1.2.inhibitor~classes$classification)$p.value,digits=3),sep=""),col = "red")
  }
}
dev.off()
