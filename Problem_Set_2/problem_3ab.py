import matplotlib.pyplot as plt

myFastq = open("I1303.1240k.fastq")

#collect read lengths from the fastq
readLengths = []
while True:
	name = myFastq.readline().strip()
	if name == "": break
	seq = myFastq.readline().strip()
	name2 = myFastq.readline().strip()
	qual = myFastq.readline().strip()
	#get the length
	curLength = len(seq)
	readLengths.append(curLength)

plt.hist(readLengths, bins = 100)

plt.savefig("readLengthHist.pdf")
		
