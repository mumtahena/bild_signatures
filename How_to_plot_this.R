# How to plot this?????

#Metadata
CCLE_Meta_data= as.matrix(read.table("~/Documents/ThesisWork/CCLE_Data/CCLE_RNASeq_Metadata.txt", sep='\t', stringsAsFactors=FALSE, header=1, row.names=1))

#determine how many cells lines for each cancer type
CCLE_MetaData_Summary = summary(CCLE_Meta_data, maxsum=31)
CCLE_MetaData_Summary
barplot(CCLE_MetaData_Summary)
CCLE_MetaData_Summary= CCLE_MetaData_Summary[1:22,1]
split_summary=do.call(rbind.data.frame, strsplit(CCLE_MetaData_Summary, ":"))
rownames(split_summary)= split_summary[,1]
colnames(split_summary)=c("Cell Line", "Freq")
split_summary=split_summary[-1]
split_summary

#barplot(split_summary[,1])
#hist(split_summary)
#as.numeric(split_summary[,1])
#hist(as.numeric(split_summary[,1]))

CCLE_cellline_frequency_file="~/Documents/ThesisWork/GitRepos/bild_signatures/CCLE_cancertype_freqs.txt"
write.table(split_summary,CCLE_cellline_frequency_file,sep='\t', col.names = NA,quote=F)
system("cat ~/Documents/ThesisWork/GitRepos/bild_signatures/CCLE_cancertype_freqs.txt ")