assign.wrapper<-function (trainingData = NULL, testData, trainingLabel, testLabel = NULL, 
          geneList = NULL, n_sigGene = NA, adaptive_B = TRUE, adaptive_S = FALSE, 
          mixture_beta = TRUE, outputDir, p_beta = 0.01, theta0 = 0.05, 
          theta1 = 0.9, iter = 2000, burn_in = 1000,control_gene=NULL) 
{
  if (is.null(geneList)) {
    pathName <- names(trainingLabel)[-1]
  }
  else {
    pathName <- names(geneList)
  }
  processed.data <- assign.preprocess(trainingData, testData, 
                                      trainingLabel, geneList, n_sigGene, theta0, theta1)
  if (!is.null(trainingData)) {
    cat("Estimating model parameters in the training dataset...\n")
    mcmc.chain.trainingData <- assign.mcmc(Y = processed.data$trainingData_sub, 
                                           Bg = processed.data$B_vector, X = processed.data$S_matrix, 
                                           Delta_prior_p = processed.data$Pi_matrix, iter = iter, 
                                           adaptive_B = FALSE, adaptive_S = FALSE, mixture_beta = TRUE)
    mcmc.pos.mean.trainingData <- assign.summary(test = mcmc.chain.trainingData, 
                                                 burn_in = burn_in, iter = iter, adaptive_B = FALSE, 
                                                 adaptive_S = FALSE, mixture_beta = TRUE)

  }
  cat("Estimating model parameters in the test dataset...\n")
  mcmc.chain.testData <- assign.mcmc(Y = processed.data$testData_sub, 
                                     Bg = processed.data$B_vector, X = processed.data$S_matrix, 
                                     Delta_prior_p = processed.data$Pi_matrix, iter = iter, 
                                     adaptive_B = adaptive_B, adaptive_S = adaptive_S, mixture_beta = mixture_beta, 
                                     p_beta = p_beta)
  mcmc.pos.mean.testData <- assign.summary(test = mcmc.chain.testData, 
                                           burn_in = burn_in, iter = iter, adaptive_B = adaptive_B, 
                                           adaptive_S = adaptive_S, mixture_beta = mixture_beta)
####Added by moom####
 pdf("Signature_convergence.pdf") 
 plot(mcmc.chain.testData$S_mcmc)
 abline(h=0,col="red")
 dev.off()
 dim=c(dim(mcmc.chain.testData$Delta_mcmc)[2],length(pathName)))
 gene_list_posterior=array(0,dim = c(dim(mcmc.chain.testData$Delta_mcmc)[2],length(pathName)))
 for(i in 1:length(pathName)){
   gene_list_posterior[,i]<-apply(mcmc.chain.testData$Delta_mcmc[,,i],2,mean)
 }
 dimnames(gene_list_posterior)=dimnames(mcmc.pos.mean.testData$Delta_pos)=dimnames(processed.data$S_matrix)
 
 deltas<-cbind(processed.data$S_matrix,processed.data$Delta_matrix,gene_list_posterior,mcmc.pos.mean.testData$Delta_pos)
 colnames(deltas)=c(paste("prior_S",pathName,sep="_"),paste("prior_delta",pathName,sep="_"),paste("posterior_mean_Delta_mcmc",pathName,sep="_"),paste("mcmc.pos.mean.testData$Delta_pos",pathName,sep="_"))
# deltas<-cbind(processed.data$S_matrix,processed.data$Delta_matrix,mcmc.pos.mean.testData$Delta_pos,mcmc.chain.testData$Delta_mcmc)
# colnames(deltas)=c("prior S","prior_delta","mcmc.pos.mean.testData$Delta_pos","mcmc.chain.testData$Delta_mcmc")

#sum(mcmc.pos.mean.testData$Delta_pos== gene_list_posterior)
 write.csv(deltas,"posterior_delta.csv",quote=F)
 #write.csv(subset(deltas,deltas[,4]>0.5,select = c(1,4)),"posterior_gene_list_delta_over_0.5.csv",quote=F)
#####End: added by moom###
  
#   if(!is.null(control_gene)){
#     plot(mcmc.chain.testData$beta_mcmc[control_gene,,])
#   }
  
  
  cat("Outputing results...\n")
  if (mixture_beta) {
    if (!is.null(trainingData)) {
      coef_train = mcmc.pos.mean.trainingData$kappa_pos
    }
    coef_test = mcmc.pos.mean.testData$kappa_pos
  }
  else {
    if (!is.null(trainingData)) {
      coef_train = mcmc.pos.mean.trainingData$beta_pos
    }
    coef_test = mcmc.pos.mean.testData$beta_pos
  }
  cwd <- getwd()
  dir.create(outputDir,showWarnings = F)##moom added this to create the output folder.
  setwd(outputDir)
  if (is.null(geneList)) {
    pathName <- names(trainingLabel)[-1]
  }
  else {
    pathName <- names(geneList)
  }
  print(paste("This analysis was run using the following parameters :",
              "n_sigGene=",n_sigGene, "adaptive_B=",adaptive_B,"adaptive_S=", adaptive_S,
              "mixture_beta=",mixture_beta,"p_beta=",p_beta,"theta0=", theta0, "theta1=",theta1, 
              "iter=",iter, "burn_in=",burn_in,"The output files are located at:",outputDir,sep=' '))###moom added this 
  if (!is.null(trainingData)) {
    rownames(coef_train) <- colnames(processed.data$trainingData_sub)
    colnames(coef_train) <- pathName
    write.csv(processed.data$S_matrix, file = "signature_gene_list_prior.csv")###moom added this to include the gene list and prior coefficient
    write.csv(coef_train, file = "pathway_activity_trainingset.csv")
  }
  rownames(coef_test) <- colnames(processed.data$testData_sub)
  colnames(coef_test) <- pathName
  write.csv(coef_test, file = "pathway_activity_testset.csv")
  if (!is.null(trainingData)) {
    heatmap.train(diffGeneList = processed.data$diffGeneList, 
                  trainingData, trainingLabel)
  }
  heatmap.test.prior(diffGeneList = processed.data$diffGeneList, 
                     testData, trainingLabel, testLabel, coef_test, geneList)
  if (adaptive_S) {
    heatmap.test.pos(testData = processed.data$testData_sub, 
                     Delta_pos = mcmc.pos.mean.testData$Delta_pos, trainingLabel, 
                     testLabel, Delta_cutoff = 0.95, coef_test, geneList)
      }
  if (!is.null(trainingData)) {
    scatter.plot.train(coef_train, trainingData, trainingLabel)
  }
  scatter.plot.test(coef_test, trainingLabel, testLabel, geneList)
  if (!is.null(testLabel)) {
    box.plot.test(coef_test, trainingLabel, testLabel, geneList)
  }
  if (!is.null(trainingData)) {
    output.data <- list(processed.data = processed.data, 
                        mcmc.pos.mean.trainingData = mcmc.pos.mean.trainingData, 
                        mcmc.pos.mean.testData = mcmc.pos.mean.testData)
  }
  else {
    output.data <- list(processed.data = processed.data, 
                        mcmc.pos.mean.testData = mcmc.pos.mean.testData)
  }
  save(output.data, file = "output.rda")
  setwd(cwd)
}
