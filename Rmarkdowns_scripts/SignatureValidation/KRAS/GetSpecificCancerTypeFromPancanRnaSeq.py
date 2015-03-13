#This script takes PANCAN20 matrix as input and outputs just the cancer type 
#you are interested in


import sys
from samUtilities import *


# Assign user input to variables/objects
RNASeqfile = file(sys.argv[1])
cancerTypeOfInterest = sys.argv[2] # if you enter "unsure" the program will print the cancer types you have
					# the options of choosing
filePathTCGA_IDconversion = sys.argv[3] # TCGA codes table report indicating cancer type for TCGA IDs
outputFile = sys.argv[4] #such as "folder/subfolder/OriginalOrSlightlyModifiedFileName"

# Get header line and determine cancer types from it
headerLine = RNASeqfile.readline().split("\t")
headerLine[-1] = headerLine[-1].rstrip("\n")

cancerTypes = getCancerTypes(headerLine[1:], filePathTCGA_IDconversion) # cancerTypes is a vector containing the cancer
									# type for each sample in the column header
									# in the same order as the samples

# If the user input "unsure" for their cancer type, print the cancer types we have
if cancerTypeOfInterest == "unsure":
        sys.exit(set(cancerTypes))

# outMatrix will contain samples with your cancer type of interest
outMatrix = []

# firstOutputLine will contain the TCGA IDs for the cancer of interest; make it as shown below
firstOutputLine = []
firstOutputLine.append(" ")

for index in range(0, len(cancerTypes)):
	if cancerTypes[index] == cancerTypeOfInterest:
		firstOutputLine.append(headerLine[index + 1])

outMatrix.append(firstOutputLine) # make the first line in the output matrix the firstOutputLine
	
# Go through each line in the input RNA-Seq file; for each line in the input
# file, add the samples 
counter = 0

for line in RNASeqfile:
	if counter % 500 == 0:
		print("Analyzing line " + str(counter))

	thisLine = line.split("\t") #get first line
	thisLine[-1] = thisLine[-1].rstrip()

	# thisOutputLine will hold this line's output
	thisOutputLine = []
	
	# Add gene symbol to the beginning of the output line
	thisOutputLine.append(thisLine[0]) 

	# Go through each item in the line and add it to the output line if it is the desired cancer type
	for index in range(1, len(thisLine)):
		if cancerTypes[index - 1] == cancerTypeOfInterest:
			thisOutputLine.append(thisLine[index])

	outMatrix.append(thisOutputLine)		
	
	counter += 1

RNASeqfile.close()

# Output outMatrix to file
outFile = file(outputFile, "w")

for line in outMatrix:
	outFile.write("\t".join(line) + "\n") 

outFile.close()
