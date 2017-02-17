import random 

#set the number of reads to simulate
numReads = 10000

#set the read length
readLength = 100

#set the genome size
genomeSize = 1000000

#initialize that we've seen each site in the genome 0 times
coveragePerSite = [0]*genomeSize

#generate coverage information for every read
for read in range(numReads):
	#draw the starting point of the read
	curStart = random.randint(0,genomeSize)
	curEnd = curStart + readLength
	#if the read is too close to the end, truncate it
	if curEnd > genomeSize: 
		curEnd = genomeSize
	#increment the number of times we've seen all the sites covered by the read
	for site in range(curStart,curEnd):
		coveragePerSite[site] += 1

#compute how many sites were seen 0 times
numSeenZeroTimes = 0
for site in coveragePerSite:
	#if the site was seen 0 times, increment the count
	if site == 0:
		numSeenZeroTimes += 1

print "The fraction of sites seen zero times is", float(numSeenZeroTimes)/genomeSize

