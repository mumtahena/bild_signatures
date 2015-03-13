# Get a mutation matrix in PANCAN20 for genes of interest

import sys
from samUtilities import *

#open files containing CNV and somatic mutation data
inputFileSomaticPath = sys.argv[1] #path for somatic mutations 1s and 0s file
genesOfInterest = sys.argv[2].split(",") #list of comma-separated genes of interest (i.e. "KRAS,PIK3CA")
filePathTCGA_IDconversion = sys.argv[3] #TCGA codes tables report indicating cancer type
					#for different XX for each TCGA-XX-YYYY
outFilePath = sys.argv[4] #COMPLETE FOLDER PATH where output files will be placed WITHOUT THE LAST FORWARD SLASH ON THE END

fileSomatic = file(inputFileSomaticPath)

#read the first line of each file to get the patient IDs
patientIDsSomatic = fileSomatic.readline().split()

print("There are " + str(len(patientIDsSomatic)) + " patients in the somatic mutation file")

#now make a somatic mutation matrix from input file for genes of interest
somaticMatrix = []
somaticFirstLine = []
somaticFirstLine.append(" ")
somaticFirstLine.extend(patientIDsSomatic)
somaticMatrix.append(somaticFirstLine)

for line in fileSomatic:
	thisLine = line.split("\t")
	thisLine[-1] = thisLine[-1].rstrip("\n")
	geneSymbol = thisLine[0]

	if geneSymbol in genesOfInterest:	
		somaticMatrix.append(thisLine)

fileSomatic.close()

# Write to file
outFile = file(outFilePath, "w")

for line in somaticMatrix:
	outFile.write("\t".join(line) + "\n")

outFile.close()
