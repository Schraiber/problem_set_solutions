import matplotlib.pyplot as plt

myFastq = open("I1303.1240k.fastq")

#initialize list of dictionaries to store base composition 
baseComposition = []
#make a dictionary of bases for every position in the reads (they are length 123 at most). Also need N
for i in range(123):
	curDict = {"A":0.0, "C":0.0, "G":0.0, "T":0.0, "N":0.0}
	baseComposition.append(curDict)

#also keep track of quality scores and how many times each position is seen
quality = [0.0]*123
timesSeen = [0.0]*123
#loop over the fastq
whichRead = 0
while True:
	if whichRead % 10000 == 0: print whichRead #this is just a simple way to track if it's running
	whichRead += 1
	name = myFastq.readline().strip()
	if name == "": break
	seq = myFastq.readline().strip()
	name2 = myFastq.readline().strip()
	qual = myFastq.readline().strip()
	#loop over the reads to get the base composition and quality
	for i in range(len(seq)):
		#get the base
		curBase = seq[i]
		#get the quality, turn it from ASCII to decimal
		curQual = ord(qual[i])-33
		#increment one additional count of curBase for this position
		baseComposition[i][curBase] += 1
		#add to the quality with the current quality score
		quality[i] += curQual
		#increment that we've seen this position one more time
		timesSeen[i] += 1

#normalize everything by dividing through by the total number of times each position is seen
for i in range(123):
	#loop over all possible nucleotides
	for nucleotide in baseComposition[i]:
		baseComposition[i][nucleotide] /= timesSeen[i]
	#also normalize quality
	quality[i] /= timesSeen[i]

#the way that I've done this, I need to extract the base composition for each site into a new list before I plot
for nuc in ["A","C","G","T","N"]:
	curComposition = []
	for i in range(123):
		curComposition.append(baseComposition[i][nuc])
	plt.plot(curComposition,label=nuc)
#put the legend in the lower right hand corner
plt.legend(loc = "lower right")
plt.savefig("baseComposition.pdf")

#make a similar figure for quality
plt.figure()
plt.plot(quality,label="quality")
plt.legend()
plt.savefig("quality.pdf")

		
