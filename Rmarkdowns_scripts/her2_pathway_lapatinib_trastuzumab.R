##HER2
subtype_her2<-pred_drug[,162]
subtype_her2=gsub("ERBB2.*","ERBB2",subtype_her2)
subtype_her2=gsub("Basal.*","BASAL",subtype_her2)
subtype_her2=gsub("Luminal.*","LUMINAL",subtype_her2)
subtype_her2=gsub(".*ormal.*","NORMAL",subtype_her2)
pred_drug[,"subtype_her2"]<-subtype_her2
her2_interesting<-subset(pred_drug,select=c(209,162,150,180))
her2_interesting[,2]<-subtype_her2

plot(her2_interesting[,3],her2_interesting[,1],main="Lapatinib response versus HER2 activity")
points(her2_interesting[her2_interesting[,2]=="NORMAL",3],her2_interesting[her2_interesting[,2]=="NORMAL",1],col="green")
points(her2_interesting[her2_interesting[,2]=="ERBB2",3],her2_interesting[her2_interesting[,2]=="ERBB2",1],col="red")
points(her2_interesting[her2_interesting[,2]=="LUMINAL",3],her2_interesting[her2_interesting[,2]=="LUMINAL",1],col="yellow")
points(her2_interesting[her2_interesting[,2]=="Claudin-low",3],her2_interesting[her2_interesting[,2]=="Claudin-low",1],col="grey")
points(her2_interesting[her2_interesting[,2]=="BASAL",3],her2_interesting[her2_interesting[,2]=="BASAL",1],col="blue")
legend("topleft",legend = c("BASAL","CLAUDIN-LOW","ERBB2","LUMINAL","NOMRAL"),pch=1,col=c("blue","grey","red","yellow","green"),cex=0.75)
trastuzumab<-read.table("~/Dropbox/Datasets/her2_amp_trastuzumab_treated_growth_rate.txt",header=1,row.names=1,sep='\t')
only_her2<-read.table("~/Dropbox/Datasets/her2_amp_lapatinib response.txt",header=1,row.names=1,sep='\t')
new_her2<-merge_drop(only_her2,data.frame(her2_interesting))
par(mar = (5,5,5,5))
plot(new_her2[,4],new_her2[,1],main="Lapatinib response versus HER2 activity\n in ERBB2 amplified breast cell lines only",pch=1,ylim=c(4.5,8),type="p",xlab="HER2 pathway status",ylab="-log(IC50 in Mols/L)",col="blue")
abline(h=5.9,col="blue",lty=3)
#points(new_her2[,4],new_her2[,1],col="dark red")
#points(her2_interesting[,3],her2_interesting[,1],pch=1,ylim=c(4.5,8),type="p",col="brown")
abline(h=5.16,col="brown",lty=3)##Lapatinib mean GI5
points(new_her2[,4],new_her2[,2],pch=1,ylim=c(4.5,7),type="p",col="brown")
abline(v=median(pred_drug[,150],na.rm=T),col="red",lty=3)
legend("topleft",legend = c("Daemen et al Lapatinib response","O'Brien et al Lapatinib response"),pch=1,col=c("blue","brown"),cex=0.7)


# plot(new_her2[,5],new_her2[,2],main="",pch=1,ylim=c(5,8),type="p",ylab = "")
# points(new_her2[,5],new_her2[,2],col="red")
# abline(h=5.92,col="red")
# legend("topleft",legend = c("Lapatinib","Trastuzumab"),pch=1,col=c("dark red","red"),cex=0.75)



plot(her2_interesting[,3],her2_interesting[,4],main="BIBW2992 response versus HER2 activity")
points(her2_interesting[her2_interesting[,2]=="NORMAL",3],her2_interesting[her2_interesting[,2]=="NORMAL",4],col="green")
points(her2_interesting[her2_interesting[,2]=="ERBB2",3],her2_interesting[her2_interesting[,2]=="ERBB2",4],col="red")
points(her2_interesting[her2_interesting[,2]=="LUMINAL",3],her2_interesting[her2_interesting[,2]=="LUMINAL",4],col="yellow")
points(her2_interesting[her2_interesting[,2]=="Claudin-low",3],her2_interesting[her2_interesting[,2]=="Claudin-low",4],col="grey")
points(her2_interesting[her2_interesting[,2]=="BASAL",3],her2_interesting[her2_interesting[,2]=="BASAL",4],col="blue")
legend("topleft",legend = c("BASAL","CLAUDIN-LOW","ERBB2","LUMINAL","NOMRAL"),pch=1,col=c("blue","grey","red","yellow","green"),cex=0.75)
abline(v=median(pred_drug[,150]),col="red",lty=2)
abline(h=6.4,col="brown",lty=2)##BIBW response threshold
par(mar=c(3.1,4.1,4.1,4.1))
boxplot(pred_drug$Sigma.AKT1.2.inhibitor~subtype_her2,ylim=c(4.5,6.5),main="Sigma AKT 1/2 inhibitor\n response in ICBP breast cancer cell lines",ylab="Sensitivity")
abline(h=median(drugs$Sigma.AKT1.2.inhibitor,na.rm=T),col="red")
par(new = T)
boxplot(pred_drug[,150]~subtype_her2,border = "grey",ylim=c(0,1), axes = F, xlab = NA, ylab = NA)
axis(side = 4)
mtext(side = 4, line = 3, "AKT pathway activity")
sigma_akt<-data.frame(Sigma.AKT1.2.inhibitor=pred_drug$Sigma.AKT1.2.inhibitor,AKT=pred_drug[,133],subtype_her2=pred_drug$subtype_her2)
#sigma_akt<-cbind(pred_drug[,150],pred_drug$Sigma.AKT1.2.inhibitor,pred_drug$subtype_her2)
#sigma_akt<-data.frame(sigma_akt[complete.cases(sigma_akt),])
#class(sigma_akt)
#colnames(sigma_akt)<-c("Sigma.AKT1.2.inhibitor", "her2","subtype_her2")
require(ggplot2)
##The values Year, Value, School_ID are
##inherited by the geoms
#dd=data.frame(Sigma.AKT1.2.inhibitor=sigma_akt[,1],her2=sigma_akt[,2],sybtype_her2=sigma_akt[,3])
png("~/akt_status.png")
ggplot(sigma_akt, aes(sigma_akt$Sigma.AKT1.2.inhibitor, sigma_akt$AKT,colour=sigma_akt$subtype_her2)) + 
  geom_point()+theme(plot.title = element_text("AKT prediction versus Sigma AKT1/2 inhibitor response"))
dev.off()

boxplot(pred_drug$GSK2141795~subtype_her2,ylim=c(0,1),main="GSK2141795 akt inhibitor\n response in ICBP breast cancer cell lines")
abline(h=median(drugs$GSK2141795,na.rm=T),col="red")
par(new = T)
boxplot(pred_drug[,150]~subtype_her2,border = "grey",ylim=c(0,1))
cor(new_her2$trastuzumab.15.ug.growth.rate.fold.change,new_her2[,4], method="spearman")
