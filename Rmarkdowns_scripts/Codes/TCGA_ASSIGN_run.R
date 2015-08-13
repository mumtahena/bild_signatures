library(data.table)
library(sva)
library(data.table)
library(devtools)
install_github("wevanjohnson/ASSIGN",ref="adapt_gene_only")
library(ASSIGN)

assign_easy_multi<-function(trainingData=train, testData=test, trainingLabel1=NULL,g=100,out_dir_base="~/Desktop/tmp",cov=0, single=0){
      if(cov==0 & single==0){
              adapB_folder<-paste(out_dir_base,paste( "adapB_multi",sep=''),sep='/')
              dir.create(file.path(out_dir_base,paste( "adapB_multi",sep='')), recursive=TRUE)
              adap_adap_folder<-paste(out_dir_base,paste( "adap_adap_multi",sep=''),sep='/')
              dir.create(file.path(out_dir_base,paste( "adap_adap_multi",sep='')), recursive=TRUE)
      }
                              
      else if (cov==0 & single==1){
              adapB_folder<-paste(out_dir_base,paste( "adapB_single",sep=''),sep='/')
              dir.create(file.path(out_dir_base,paste( "adapB_single",sep='')), recursive=TRUE)
              adap_adap_folder<-paste(out_dir_base,paste( "adap_adap_single",sep=''),sep='/')
              dir.create(file.path(out_dir_base,paste( "adap_adap_single",sep='')), recursive=TRUE)
     }
                                                        
                                                          
     set.seed(1234)
     assign.wrapper(trainingData=trainingData, testData=testData, trainingLabel=trainingLabel1, geneList=NULL, n_sigGene=g, adaptive_B=T, adaptive_S=F, mixture_beta=F, outputDir=adapB_folder, theta0=0.05, theta1=0.9, iter=100000, burn_in=5000)  
                                                                
     set.seed(1234)
     assign.wrapper(trainingData=trainingData, testData=testData, trainingLabel=trainingLabel1, geneList=NULL, n_sigGene=g, adaptive_B=T, adaptive_S=T, mixture_beta=F, outputDir=adap_adap_folder, theta0=0.05, theta1=0.9, iter=100000, burn_in=5000) 
                                                                      
}
