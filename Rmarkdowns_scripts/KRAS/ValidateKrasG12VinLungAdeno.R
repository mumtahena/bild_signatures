# Takes as input a folder containing folders that have
# probabilities files from BinReg runs analyzing TCGA lung adeno samples; also takes as input
# a file containing K-Ras mutant status for for those sample lung adeno samples.

# From this input the script generates a p-value (in the command-line) telling you 
# how well the model discriminated between K-Ras mutant and WT lung adeno samples (simple
# t-test of probabilities between WT and mutant samples). Further, the script generates
# output that you can use to plot the differences.

# Get command-line arguments - 1st argument is directory with subdirectories containing
# BinReg probabilities.txt files, and 2nd argument is a file containing mutation
# status for each patient for the gene of interest (patients in columns on first row; 
# second row contains 1s or 0s, 0 for WT and 1 for mutant).
args <- commandArgs(TRUE)
directoryOfBinRegOutputDirectories <- args[1]
mutationStatusFile <- args[2]
outputPDFfile <- args[3] # user must also input the output pdf file 

# Output to PDF
pdf(outputPDFfile)

# Get original directory name
originalDirectory <- getwd() 

# Read in mutation status information from mutationStatusFile
mutationStatus <- as.data.frame(read.table(mutationStatusFile, header=TRUE, stringsAsFactors=FALSE, row.names=1, check.names=FALSE))

# Go to directory with BinReg output and get all the subdirectory names
setwd(directoryOfBinRegOutputDirectories)
foldersWithProbabilities <- list.files() # get subdirectory names (no files should be in this folder, only subdirectories)

cat("\nI will look for BinReg probabilities files in the following folders: ")
cat(foldersWithProbabilities)
cat("\n")

# Go through each probabilities file in each subdirectory and
# compare those K-Ras probabilities with the K-Ras mutation status
# and output the data as a t-test to the command line and as a graph
# in a .pdf
for (folder in foldersWithProbabilities)
{
	# Go to subfolder
	setwd(folder)

	# Let the user know what's going on
	cat("\nNow in ")
	cat(getwd())
	cat("\n\tLooking at the probabilities.txt file")


	# Read in the probabilities.txt file into a dataframe
	probabilities <- as.data.frame(read.table("probabilities.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t"))
	
	# These variables will be vectors containing the probabilities in WT and K-Ras mutant samples, respectively
	probabilitiesKrasWt <- NULL
	probabilitiesKrasMut <- NULL

	# Go through the probabilities.txt row-by-row and add probabilities to the appropriate vector defined above
	numRows <- nrow(probabilities) 	

	for (rowNum in 1:numRows)
	{
		# First make sure the line starts with "TCGA" (i.e. is not a training sample) before considering it
		if (substr(probabilities[rowNum, 2], 1, 4) == "TCGA")
		{
			# Get the first 12 characters of the TCGA ID
			TcgaIdMinimal <- substr(probabilities[rowNum, 2], 1, 12)
	
			# Find the column index where the patient is found in the mutation data.frame
			columnIdThisPatient <- grep(TcgaIdMinimal, colnames(mutationStatus))
			
			# Now if the columnId was found (i.e. if its length is not zero) add the probability for the patient
			# that column to the appropriate vector (WT or mutant) based on the mutation status (1 or 0) in the columnId
			if (length(columnIdThisPatient) != 0)
			{
				# Find out if this patient is Kras mutant (1) or WT (0)
				isThisPatientKrasMutant <- mutationStatus[1, columnIdThisPatient]

				# If mutant, add the probability to the probabilitiesKrasMut vector, otherwise add to the WT vector
				if (isThisPatientKrasMutant == 1)
				{
					# probabilities are in column 5, just FYI
					probabilitiesKrasMut <- c(probabilitiesKrasMut, probabilities[rowNum, 5])
				}
				else
				{
					probabilitiesKrasWt <- c(probabilitiesKrasWt, probabilities[rowNum, 5])
				}
			}
		}
	}

	numWT <- length(probabilitiesKrasWt)
	numMutant <- length(probabilitiesKrasMut)

	cat("\n\tMean of WT probabilities is ")
	cat(mean(probabilitiesKrasWt))
	
	cat("\n\tMean of mutant probabilities is ")
	cat(mean(probabilitiesKrasMut))

	cat("\n\tT-test for this pair is ")
	pVal <- t.test(probabilitiesKrasWt, probabilitiesKrasMut)$p.value
	cat(pVal)
	cat("\n")

	parameters <- strsplit(folder, "_")

	stripchart(list(probabilitiesKrasWt, probabilitiesKrasMut), vertical = TRUE, method = "jitter", pch = 21, col = "maroon", bg = "bisque", cex.lab=1, cex.axis=1, xaxt="n", xlab="Lung adenocarcinoma K-Ras mutation status", ylab="Probability of K-Ras (G12V) activation", xlim=c(0,3), ylim=c(0,1))

	axis(1, at=1:2, labels=c(paste0("WT (n=", numWT, ")"), paste0("Mutant (n=", numMutant, ")")))
	
	boxplot(probabilitiesKrasWt, probabilitiesKrasMut, add=TRUE, axes=FALSE, outline=FALSE, xlim=c(0,3), ylim=c(0,1))
	title(paste0("K-Ras G12V signature validation\nParameters: ", parameters[[1]][1], " normalization, ", parameters[[1]][3], " genes, ", parameters[[1]][2], " metagene(s)", "\n(p=", sprintf("%.3f", pVal), ")"), cex=0.8)

	# Go up one directory so you can go into the next subfolder	
	setwd("../")	
}

dev.off()
