#Utilities file for functions

def getCancerTypes(IDlist, IDconversionFile):
	#make a dictionary based on the TCGA ID conversion file
	TCGAdictionary = {}
	for line in file(IDconversionFile):
		thisLine = line.rstrip().split("\t")
		TCGAdictionary[thisLine[0]] = thisLine[1]
		
	#print(TCGAdictionary)
	#print("DID IT PRINT?\n")
	
	cancerTypeList = []
	
	for item in IDlist:
		if ("TCGA-" in item):
			cancerTypeList.append(TCGAdictionary[item[5:7]])
	
	#print(set(cancerTypeList))
	
	return cancerTypeList

def isNumber(s):
	try:
		float(s)
		return True
	except ValueError:
		return False