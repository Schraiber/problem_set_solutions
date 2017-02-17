def read_fasta(fastaFile):
	#initialize empty dictionary
	fastaDict = {}
	#go over every line int he file
	for line in fastaFile:
		#if it's got a > at the start, we know it's the name of a sequence
		if line[0] == ">":
			chromName = line.strip()[1:]
			#initialize newly seen sequence names with empty lists 
			#this way, we can append every line as it comes in
			fastaDict[chromName] = []
		else:
			#append all the lines that belong to that chromosome
			fastaDict[chromName].append(line.strip())
	#now we need to merge all the lines together
	#loop over every chromosome and merge all of its lines
	for chrom in fastaDict:
		#''.join(list_of_strings) will join all the strings in the list into a new string with no separating characters
		fastaDict[chrom] = ''.join(fastaDict[chrom])
	return fastaDict

#helper function to give complementary nucleotides
def complement(nuc):
	if nuc == "A":
		return "T"
	elif nuc == "C":
		return "G"
	elif nuc == "G":
		return "C"
	elif nuc == "T":
		return "A"

def reverse_complement(seq):
	#first reverse the sequence using a weird python trick
	revSeq = seq[::-1]
	#now go through and complement everything
	revCom = []
	for pos in revSeq:
		curComplement = complement(pos)
		revCom.append(curComplement)
	#now join the reverse complement back together
	revComSeq = ''.join(revCom)
	return revComSeq

def naiveStringSearch(substring,string):
	#loop over all positions in the genome (but don't need to consider positions where you would have the substring fall off the edge)
	for i in range(len(string)-len(substring)):
		#go through all positions in the substring
		for j in range(len(substring)):
			#if we find a place where it mismatches, we know we don't have a match, so break!
			if substring[j] != string[i+j]:
				break
			#if we get to the end of the substring and haven't broken
			#we must have found a match!
			if j == len(substring)-1:
				return i #NB: a return IMMEDAITELY exits the function
	#reutnr None if it failed
	return None

myGenome = open("my_genome.fa")

myGenomeDict = read_fasta(myGenome)

#reverse complmement the genome
myGenomeComplementDict = {}
for chrom in myGenomeDict:
	myGenomeComplementDict[chrom] = reverse_complement(myGenomeDict[chrom])

myFastq = open("my_reads.fastq")

myOutput = open("string_matching_mapping.txt","w")

whichRead = 0
while True:
	if whichRead % 10000 == 0: print whichRead
	whichRead += 1
	#read an entry the normal way
	name = myFastq.readline().strip()
	if name == "": break
	seq = myFastq.readline().strip()
	name2 = myFastq.readline().strip()
	qual = myFastq.readline().strip()
	#go over all the chromosomes
	for chrom in myGenomeDict:
		#first search forward chromosome
		mapping = naiveStringSearch(seq,myGenomeDict[chrom])
		if mapping is not None:
			#this means we found it!
			myOutput.write("%s\t%d\t+\t%s\n"%(chrom,mapping,name))
			#we can break because we don't need to search any more
			break 
		#then search reverse
		mapping = naiveStringSearch(seq,myGenomeDict[chrom])
		if mapping is not None:
			#this means we found it!
			myOutput.write("%s\t%d\t-\t%s\n"%(chrom,mapping,name))
			#we can break because we don't need to search anymore
			break
