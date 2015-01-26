# ---
# title: "IGF1R"
# author: "Shelley"
# date: "December 16, 2014"
# output: html_document
# ---


```{r}
library(sva)
library(ASSIGN)
```


Read in the clinical data

preTreatment<- as.matrix(read.table("~/Downloads/TCGA_Breast_PreTreatmentStatus_NA_Removed.txt", row.names=1, header=FALSE, check.names=FALSE, stringsAsFactors=F))
View(preTreatment)
dim(preTreatment)

Read in the datasets
```{r}
setwd("~/Dropbox/bild_signatures/Datasets")

GFP18_AKT_BAD_HER2_IGF1R_RAF_ERK <- as.matrix(read.table("Documents/ThesisWork/GitRepos/bild_signature_validation_old_repo/Datasets/GFP18_AKT_BAD_HER2_IGF1R_RAF_ERK.tpmlog", stringsAsFactors=FALSE, row.names=1))

head(GFP18_AKT_BAD_HER2_IGF1R_RAF_ERK)
row.names(GFP18_AKT_BAD_HER2_IGF1R_RAF_ERK)
colnames(GFP18_AKT_BAD_HER2_IGF1R_RAF_ERK)

TCGA_Breast_RNAseq <- read.delim("~/Documents/ThesisWork/GSOA_Manuscript/GSOA_Files/RNA_Seq_Files/PANCAN12.IlluminaHiSeq_RNASeqV2.geneExp.tumor_whitelist_Breast_nodup", row.names=1, check.names=FALSE, stringsAsFactors=F)
TCGA_Breast_RNAseq=as.matrix(TCGA_Breast_RNAseq)

head(TCGA_Breast_RNAseq)
length(colnames(TCGA_Breast_RNAseq))
colnames(TCGA_Breast_RNAseq)
row.names(TCGA_Breast_RNAseq)

gfp_BAD <-subset(GFP18_AKT_BAD_HER2_IGF1R_RAF_ERK,select=c(GFP.1:GFP.12,BAD.1:BAD.6))
head(gfp_BAD)
length(colnames(gfp_IGF1R))

head(gfp_BAD[,1:17])

# the 1 means by row name (gene names), if ~80% are zero then we don't want that gene 
gfp_BAD <- gfp_BAD[apply(gfp_BAD[,1:17]==0,1,mean) < 0.8,]###if not removed over 83%; error occurs in batch adjustment
nrow(gfp_BAD)
removed= nrow(gfp_BAD)-nrow(gfp_BAD)
print( "number of genes removed after filtering = ") 
removed
```

Some function to make life easier :)
```{r}
merge_drop<-function(x,y,by=0)
{
  new_m<-merge(x,y,by=by)
  rownames(new_m)<-new_m$Row.names
  return(new_m[,2:length(colnames(new_m))])
}

pcaplot<-function(mat,sub,scale=T){
  if(sum(sub)!=length(mat)){
    print("verify the subscripts...exiting now")
    }
  else{
    pca_mat <- prcomp(t(mat), center=T,scale=scale)
    plot(pca_mat)
    plot(pca_mat$x[,1],pca_mat$x[,2])
    index= 1
    for(i in 1:length(sub)){
      #print(rownames(pca_mat$x)[index:sub[i]+index-1],"has color", )
      print(i)
      if(i==1){
        points(pca_mat$x[1:sub[i]],pca_mat$x[1:sub[i],2],col=i+1)
        }
      else if(i==length(sub)){
         points(pca_mat$x[index:length(rownames(pca_mat$x))],pca_mat$x[index:length(rownames(pca_mat$x)),2],col=i+1)
         }
      else{
        points(pca_mat$x[index:index+sub[i]],pca_mat$x[index:index+sub[i],2],col=i+1)
        }
       index=index+sub[i]
      }
  }
}
```


# Can use the commented out ones or not. Filtering out lowered expression or low varything genes 
# With adaptive and non-adaptive you don't really need too. It might make a different with non-adaptive 
```{r, cache=TRUE}
#AVG_IN <- apply(gfp_her2_1[,1:17], 1, mean)
#VAR_IN <- apply(gfp_her2_1[,1:17], 1, var)
#threshold_E <- sort(AVG_IN)[round(length(AVG_IN)*0.1)]
#threshold_V <- sort(VAR_IN)[round(length(VAR_IN)*0.1)]
#gfp_her2_1<- gfp_her2_1[which((AVG_IN > threshold_E)&(VAR_IN > threshold_V)),]
#nrow(gfp_her2_1)
# merging the test and traning data in order to use combat, need to be in the same matrix

dim(gfp_BAD) #177795 18
head(gfp_BAD)
dim(TCGA_Breast_RNAseq) # 508
head(TCGA_Breast_RNAseq)
gfp_BAD_breast<-merge_drop(gfp_BAD,TCGA_Breast_RNAseq ,by=0)
dim(gfp_BAD_breast) # 16165 526
head(gfp_BAD_breast)
# specify the number of gfp,bad, and breast
sub<-c(12,6,508)
#creates PCA to look at batch differences
pcaplot(mat =gfp_BAD_breast,sub=sub)
```


Running combat without covariates
```{r}
##############Batch adjustment is needed############
# creating a matrix that tells combact which ones comes from which batch
bat2<-as.matrix(cbind(c(colnames(gfp_BAD),colnames(TCGA_Breast_RNAseq)),c(rep(1,length(colnames(gfp_BAD))),rep(2,length(colnames(TCGA_Breast_RNAseq))))))
bat2
bat2[,2]
# calls combat, and batch says what the batches for each sample are, usually always do "Null", using diff covariates did not get better results 
combat_expr2<-ComBat(dat=gfp_BAD_breast, batch=bat2[,2], mod=NULL, numCovs=NULL)
head(combat_expr2)
# want to see them mixed together
pcaplot(combat_expr2,sub)
```

Run ASSIGN with no covariates and different number of genes
```{r}
##########calling ASSIGN########
###making a easy to call ASSIGN function#######
assign_easy<-function(trainingData=train, testData=test, trainingLabel1=NULL,g=150,out_dir_base="~/Desktop/tmp/",cov=0){
if(cov==0){
  adap_folder<-paste(out_dir_base,paste( "adap",g,sep=''),sep='/')
  dir.create(file.path(out_dir_base,paste( "adap",g,sep='')))
  nonadap_folder<-paste(out_dir_base,paste( "nonadap",g,sep=''),sep='/')
  dir.create(file.path(out_dir_base,paste( "nonadap",g,sep='')))
  }
else{
  adap_folder<-paste(out_dir_base,paste( "adap_cov",g,sep=''),sep='/')
  dir.create(file.path(out_dir_base,paste( "adap_cov",g,sep='')))
  nonadap_folder<-paste(out_dir_base,paste( "nonadap_cov",g,sep=''),sep='/')
  dir.create(file.path(out_dir_base,paste( "nonadap_cov",g,sep='')))
}

set.seed(1234)
assign.wrapper(trainingData=trainingData, testData=testData, trainingLabel=trainingLabel1, testLabel=NULL, geneList=NULL, n_sigGene=g, adaptive_B=T, adaptive_S=F, mixture_beta=F, outputDir=adap_folder, theta0=0.05, theta1=0.9, iter=2000, burn_in=1000)  

set.seed(1234)
assign.wrapper(trainingData=trainingData, testData=testData, trainingLabel=trainingLabel1, testLabel=NULL, geneList=NULL, n_sigGene=g, adaptive_B=F, adaptive_S=F, mixture_beta=F, outputDir=nonadap_folder, theta0=0.05, theta1=0.9, iter=2000, burn_in=1000)  

}

ncol(combat_expr2) #526
trainingLabel1 <- list(control = list(bad=1:12), bad=13:18) 
trainingLabel1
# what are the labels for the training set
train_bad<-combat_expr2[,1:18] # 18 is training
head(train_bad)
test<-combat_expr2[,19:526] # the test
head(test)
ncol(test)
ncol(train_bad)
dim(test)
head(combat_expr2[19])


#create the output directory for ASSIGN 
dir.create("~/Desktop/tmp/bad")
setwd("~/Desktop/tmp/bad")
assign_easy(trainingData=train_bad, testData = test,trainingLabel1=trainingLabel1,g=100,out_dir_base="~/Desktop/tmp/bad", cov=0)



```

#read back in the adaptive results
bad_predctions_adaptive<-read.csv(file="adap100/pathway_activity_testset.csv", header=1,row.names=1) 
View(bad_predctions_adaptive)
summary(bad_predctions_adaptive)

bad_predctions_nonadaptive<-read.csv(file="nonadap100/pathway_activity_testset.csv", header=1,row.names=1) 
View(bad_predctions_nonadaptive)
summary(bad_predctions_nonadaptive)
class(bad_predctions_nonadaptive)
cor(bad_predctions_adaptive[1],bad_predctions_nonadaptive[1])
# Split the predictions in high and low

bad_predctions_nonadaptive[1]

#sub set based on high and low'


# adaptive
bad_on_adp_mean=subset(bad_predctions_adaptive, bad>0.463)
dim(bad_on_adp_mean)

bad_low_adp_mean=subset(bad_predctions_adaptive, bad<0.463)
dim(bad_high_adp_mean)

bad_off_adp_5=subset(bad_predctions_adaptive, bad<0.5)
dim(bad_off_adp_5)

bad_on_adp_5=subset(bad_predctions_adaptive, bad>0.5)
dim(bad_on_adp_5)


# nonadaptive
bad_on_nonadp_mean=subset(bad_predctions_nonadaptive, bad>0.38)
dim(bad_on_nonadp_mean)

bad_off_nonadp_mean=subset(bad_predctions_nonadaptive, bad<0.38)
dim(bad_on_nonadp_mean)

bad_off_nonadp_5=subset(bad_predctions_nonadaptive, bad<0.5)
head(bad_off_nonadp_5) #440 with bad low
dim(bad_off_nonadp_5)

bad_on_nonadp_5=subset(bad_predctions_nonadaptive, bad>0.5)
head(bad_on_nonadp_5) #68 with bad high 
dim(bad_on_nonadp_5)

#since most of these are pretreatment samples, bad seems to be off more commonly than on. 

#merge preTreament with predictions

head(preTreatment)
predicts_pretreat_nonadap=merge_drop(bad_predctions_nonadaptive, preTreatment, by=0)
View(predicts_pretreat_nonadap)
YES=subset(predicts_pretreat_nonadap, predicts_pretreat_nonadap[2]=="YES")
head(YES)

NO=subset(predicts_pretreat_nonadap, predicts_pretreat_nonadap[2]=="NO")
NO



# EMT diff expression anaysis

EMTgenes = TCGA_Breast_RNAseq[c("VIM", "SNAI1", "TWIST1", "TGFB1", "CDH2", "ZEB1"),]
head(EMTgenes)
dim(EMTgenes)
# Get the mean of the 6 genes
Gene_expr_sum= rowSums(EMTgenes,na.rm = FALSE, dims = 1)
View(Gene_expr_sum)
row.names(Gene_expr_sum)
sum_all_6=mean(Gene_expr_sum[1])
View(sum_all_6)

mean(Gene_expr_sum[1], na.rm = TRUE)
sum=(10784813.54/508)+(67181.97/508)+(68318.85/508)+(509022.12/508)+(143736.54/508)+(365548.61/508)
sum

mean=sum/6
mean

#mean of the total expression
mean_all_6genes=sum/508
mean

#mean of the sample expression
mean_samples=mean_all_6genes/508
mean_samples #3916

#get the mean on the 6 genes in the samples
samples_expr_sum= colMeans(EMTgenes,na.rm = FALSE, dims = 1)
samples_expr_sum=as.data.frame(samples_expr_sum)
class(samples_expr_sum)
head(samples_expr_sum[1])
View(samples_expr_sum)
class(samples_expr_sum)
samples_expr_sum[1]

summary(samples_expr_sum)

EMT_high=subset(samples_expr_sum, samples_expr_sum>3398.9)
dim(EMT_high)

EMT_low=subset(samples_expr_sum, samples_expr_sum>3398.9)
dim(EMT_low)

emtmatrix<-as.matrix(cbind(c(row.names(EMT_high),row.names(EMT_low)),c(rep(1,length(row.names(EMT_high))),rep(2,length(row.names(EMT_low))))))
emtmatrix[,2]


# run the diff exp anaylsis.
library(limma)
library(GSVA)
library(gage)

design <- model.matrix(~ factor(emtmatrix[,2]))
View(design)
colnames(design) <- c("ALL", "EMTvsN")

fit <- lmFit(TCGA_Breast_RNAseq, design)
head(fit)
fit <- eBayes(fit)
View(fit)

EMTgenes_df=topTable(fit, coef="EMTvsN", number=Inf)
View(EMTgenes)
write.table(EMTgenes_df, file = "DF_expr_analysis_EMThighvsLow.txt", sep = "\t", col.names=NA)
