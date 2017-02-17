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

def create_string_hash(string,substring_length):
	#initialize empty dictionary
	stringHash = {}
	#go over everything up until substrings would be falling off the edge
	for i in range(len(string)-substring_length):
		#extract the substring that goes from position i to position i+substring_length
		curSubstring = string[i:(i+substring_length)]
		#add it as a key to the hash, where the value is the position
		stringHash[curSubstring] = i
	return stringHash

def search_hash(substring, stringHash):
	if substring in stringHash:
		#if it's in there, just return the value, which is the position
		return stringHash[substring]
	else:
		return None
	

myGenome = open("my_genome.fa")

myGenomeDict = read_fasta(myGenome)

#reverse complmement the genome
myGenomeComplementDict = {}
for chrom in myGenomeDict:
	myGenomeComplementDict[chrom] = reverse_complement(myGenomeDict[chrom])

#preprocess genome hash
myGenomeHashed = {}
for chrom in myGenomeDict:
	myGenomeHashed[chrom] = create_string_hash(myGenomeDict[chrom],100)

#preprocess reverse genome hash
myGenomeComplementHashed = {}
for chrom in myGenomeComplementDict:
	myGenomeComplementHashed[chrom] = create_string_hash(myGenomeComplementDict[chrom],100)

myFastq = open("my_reads.fastq")

myOutput = open("hashing_mapping.txt","w")

whichRead = 0
while True:
	if whichRead % 500 == 0: print whichRead
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
		mapping = search_hash(seq,myGenomeHashed[chrom])
		if mapping is not None:
			#this means we found it!
			myOutput.write("%s\t%d\t+\t%s\n"%(chrom,mapping,name))
			#we can break because we don't need to search any more
			break 
		#then search reverse
		mapping = search_hash(seq,myGenomeComplementHashed[chrom])
		if mapping is not None:
			#this means we found it!
			myOutput.write("%s\t%d\t-\t%s\n"%(chrom,mapping,name))
			#we can break because we don't need to search anymore
			break
