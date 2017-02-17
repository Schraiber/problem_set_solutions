def parse_fasta(fastaFile):
	#initialize an empty dictionary
	fastaDict = {}
	for line in fastaFile:
		if line[0] == ">":
			#if the first charaacter is >, then we know it's a sequence name
			contigName = line.strip()
		else:
			#otherwise it's the sequence, so we want to assign it to that entry in the dictionary
			contigSequence = line.strip()
			fastaDict[contigName] = contigSequence
	return fastaDict

def computeNx(contigs,x):
	#get the lengths of all the contigs
	contigLengths = []
	for contig in contigs:
		contigLengths.append(len(contigs[contig]))
	#sort the list of contig lengths
	sortedLengths = sorted(contigLengths)
	#generate a new list that has repeats of all the sorted lengths, as outlined in the lecture notes
	repeatedList = []
	for length in sortedLengths:
		repeatedList.extend([length]*length)
	#get the 1-x quantile of the list
	pos = int((1.0-x)*len(repeatedList))
	return repeatedList[pos]

myFastaFile = open("test.fasta")

myFastaDict = parse_fasta(myFastaFile)

#generate output and write it to a file
#NB: This is hilariously inefficient, since I recreate the same exact lists every single time
outFile = open("problem_2_output.txt", "w")
x = 0.05
while x < 1:
	Nx = computeNx(myFastaDict,x)
	outFile.write("x = %f, Nx = %d\n"%(x,Nx))
	x += 0.05
	
