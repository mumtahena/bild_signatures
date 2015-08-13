if (!require("plyr")) {
  install.packages("plyr", dependencies = TRUE)
  library(plyr)
}

source('~/Dropbox/bild_signatures/bild_signatures/Rmarkdowns_scripts//Key_ASSIGN_functions.Rmd', echo=F)
setwd("~/Dropbox/gdsc_drugs/")
filenames<-system("ls *.csv", intern=TRUE)
filenames
data_gdsc=NULL
for(i in 1:length(filenames))
{
  f<-read.csv(filenames[i], header=1) ###reading in the filess one at a time
  rownames(f)<-f[,2]
  f_select<-subset(subset(f,f[,4]=="breast"),select=c(7))
  colnames(f_select)<-levels(f[1,1])[1]
  print(filenames[i])
  print( colnames(f_select))
  if(i==1){
    data_gdsc<- f_select
    print(dim(data_gdsc))
  }
  else{
    data_gdsc<- merge_drop(data_gdsc,f_select,by=0,all = T)
    print(dim(data_gdsc))
  }
}
head(data_gdsc)
rownames(data_gdsc)=gsub("-","",rownames(data_gdsc))
colnames(data_gdsc)=gsub("-","",colnames(data_gdsc))
data_gdsc<-data_gdsc[,order(colnames(data_gdsc))]
head(data_gdsc)
write.table(data_gdsccor(icbp_gdsc$Lapatinib.x,icbp_gdsc$Lapatinib.y,use="pairwise")
,"GDSC_drugs_breast_cancer_cell_lines.txt",sep='\t',col.names=NA,quote=F)
colnames(data_gdsc)
