# This scripts analyzes the difference bettwen MCL-1 dependent and independent cell lines 

setwd("~/Documents/Multipathway_Modeling/")
source("http://www.r-statistics.com/wp-content/uploads/2010/02/Barnard.R.txt")
source("http://www.r-statistics.com/wp-content/uploads/2012/01/source_https.r.txt")
source_https("https://raw.github.com/talgalili/R-code-snippets/master/Barnard.R")
MCLdata=read.csv("Mcl1_Xiao_Data.txt",header = T,sep='\t')[1:23,1:5]

View(MCLdata)
colnames(MCLdata)

MCLdata[1:13,4]


boxplot(MCLdata[1:14,4],MCLdata[15:23,4], ylab="% Viability upon treatment (MCL-1 siRNA)",names = c("AKT Phenotype \n n=14", " EGFR Phenotype \n n=9"), col = c("red", "grey"), main="Mcl-1 Dependence Across Growth Phenotypes", ylim=c(0,100))
text(1,2,"p-value = 0.095")
t.test(MCLdata[1:14,4],MCLdata[15:23,4])


MclDependence<- matrix(c(9,2,5,7), nrow=2, dimnames=list(c("EGFR", "AKT"), c("Mcl Dependent", "Non-Mcl Dependent")))
fisher.test(MclDependence)
fisher.test(MclDependence, alternative = "greater")
View(MclDependence)
Barnard(MclDependence)


TeaTasting <-
  matrix(c(3, 1, 1, 3),
         nrow = 2,
         dimnames = list(Guess = c("Milk", "Tea"),
                         Truth = c("Milk", "Tea")))
View(TeaTasting)
fisher.test(TeaTasting, alternative = "greater")

Convictions <-
  matrix(c(2, 10, 15, 3),
         nrow = 2,
         dimnames =
           list(c("Dizygotic", "Monozygotic"),
                c("Convicted", "Not convicted")))
View(Convictions)



