# This script runs a fisher test on the quanified western blots results for the 
# pathway modeling manuscript. 

pBad<- matrix(c(9,2,5,7), nrow=2, dimnames=list(c("EGFR", "AKT"), c("pBad Expressed", "pBad Not-Expressed")))

Bim_light<- matrix(c(4,7,6,5), nrow=2, dimnames=list(c("EGFR", "AKT"), c(("BIM_Light Expressed", "BIM_Light Not-Expressed")))
Bim_dark<- matrix(c(4,6,11,8), nrow=2, dimnames=list(c("EGFR", "AKT"), c("BIM_Dark Expressed", "BIM_Dark Not-Expressed")))

fisher.test(pBad)
fisher.test(Bim_light)
fisher.test(Bim_dark)