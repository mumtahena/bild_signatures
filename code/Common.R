plotSimulatedResults = function(numSignal, metricValues, outFilePath, ylab)
{
  data = as.data.frame(cbind(numSignal, metricValues))
  data = data[order(data$metricValues),]

  pdf(outFilePath, height=4, width=6)
  boxplot(data$metricValues ~ as.factor(data$numSignal), ylim=c(0, 1), xlab="# Signal Genes", ylab=ylab, main="")
  graphics.off()
}

plotSimulatedAUCByNumGenes = function(numSignal, metricValues, outFilePath, ylab)
{
  data = as.data.frame(cbind(numSignal, metricValues))

  numSignalOptions = sort(unique(numSignal))
  numSignalOptions = numSignalOptions[which(numSignalOptions>0)]

  aucs = NULL
  for (numSignalOption in numSignalOptions)
  {
    data2 = data[which(data$numSignal %in% c(numSignalOption, 0)),]
    aucs = c(aucs, calcAUCIntervals(data2$metricValues, data2$numSignal > 0, direction=">"))
    print(numSignalOption)
    print(aucs)
  }
  stop()

  pdf(outFilePath, height=4, width=6)
  boxplot(data$metricValues ~ as.factor(data$numSignal), ylim=c(0, 1), xlab="# Signal Genes", ylab=ylab, main="")
  graphics.off()
}

calcAUC = function(numericValues, answers, direction="<")
{
  library(pROC)
  roc_result = roc(formula = as.factor(answers) ~ numericValues, ci=TRUE, plot=FALSE, print.auc=FALSE, main=markers, direction=direction)

  return(roc_result$ci[2])
}

plotROC = function(outFilePath, classValues, predictionValues)
{
  library(pROC)
  pdf(outFilePath)
  roc_result = roc(formula = as.factor(classValues) ~ predictionValues, ci=TRUE, plot=FALSE, print.auc=FALSE, main="", direction="auto")
  plot(1 - roc_result$specificities, roc_result$sensitivities, type="l", xlab="1 - Specificity", ylab="Sensitivity", lwd=3)
  abline(a=0, b=1, col="gray")
  auc = roc_result$ci[2]
  text(0.4, 0.4, paste("AUC: ", auc, sep=""))
  graphics.off()
}
