
#ASSIGN single pathway
assign_easy<-function(trainingData=train, testData=test, trainingLabel1=NULL,g=100,out_dir_base="~/Desktop/tmp/",cov=0){
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
  
  #set.seed(1234)
  #assign.wrapper(trainingData=trainingData, testData=testData, trainingLabel=trainingLabel1, testLabel=NULL, geneList=NULL, n_sigGene=g, adaptive_B=T, adaptive_S=F, mixture_beta=F, outputDir=adap_folder, theta0=0.05, theta1=0.9, iter=2000, burn_in=1000)  
  
  set.seed(1234)
  assign.wrapper(trainingData=trainingData, testData=testData, trainingLabel=trainingLabel1, testLabel=NULL, geneList=NULL, n_sigGene=g, adaptive_B=F, adaptive_S=F, mixture_beta=F, outputDir=nonadap_folder, theta0=0.05, theta1=0.9, iter=10000, burn_in=1000)  
}
