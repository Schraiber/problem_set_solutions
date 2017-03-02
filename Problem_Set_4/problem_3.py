import pysam
import matplotlib.pyplot as plt

#open the bam file
bamFile = pysam.AlignmentFile("BC187_mapped.sorted.bam")

#initialize some empty lists to store all the things
mappingQualities = []
strands = []
for read in bamFile.fetch():
	#get the mapping quality of the read
	mappingQuality = read.mapping_quality
	#get the strand that the read maps to
	#True means - strand, False means + strand
	strand = read.is_reverse
	mappingQualities.append(mappingQuality)
	strands.append(strand)

#histogram the mapping qualities
plt.hist(mappingQualities)
#save it
plt.savefig("mapping_quality.pdf")

#Compute the proportion of reads that are + strand
number_on_plus_strand = 0.0
for strand in strands:
	if strand == False: number_on_plus_strand += 1
frac_on_plus = number_on_plus_strand/len(strands)
print "%f of the reads are on the + strand"%(frac_on_plus)
