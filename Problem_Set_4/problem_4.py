import pysam
from scipy.special import binom

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

#write a separate function to compute genotype likelihoods
def compute_genotype_likelihood(n, alt_num,e=0.01):
	#computing n will make the code cleaner
	n = alt_num + ref_num 
	#use binom from scipy.special. I googled it.
	GLref = binom(n,alt_num)*e**alt_num*(1.0-e)**ref_num
	GLalt = binom(n,alt_num)*(1.0-e)**alt_num*e**ref_num
	return GLref, GLalt

bamFile = pysam.AlignmentFile("BC187_mapped.sorted.bam")

#open and parse the reference genome
fastaFile = open("yeast_genome.fa")
refGenome = read_fasta(fastaFile)

outFile = open("yeast_genotype_likelihoods.txt","w")

#loop over all positions with coverage
for pileupcolumn in bamFile.pileup():
	#need to get chromosome and position to check the reference genome
	chrom = pileupcolumn.reference_name
	pos = pileupcolumn.pos
	#get the reference genome
	ref_allele = refGenome[chrom][pos]
	#get the total coverage at this site
	n = pileupcolumn.n
	#skip the site if coverage is too high
	if n > 100: continue
	#going to need to get the number of reference alleles seen
	ref_num = 0
	#we could conceivably see all possible alleles, so a dictionary is a good data structure
	alt_alleles = {}
	#loop over the reads at that position
	for pileupread in pileupcolumn.pileups:
		#skip sites that are indels
		if pileupread.is_del or pileupread.is_refskip: continue
		#get the nucleotide from that read
		read = pileupread.alignment.query_sequence[pileupread.query_position]
		#if it's not the reference allele, need to add to dictionary or incrememnt
		if read != ref_allele:
			if read in alt_alleles:
				#if it's in the dictionary, increment
				alt_alleles[read] += 1
			else:
				#if it's not, add one
				alt_alleles[read] = 1
		else:
			#the allele is ref, so count it as one more time seeing ref
			ref_num += 1
	#go through and get the alternate allele that's seen the most times
	#before I've seen any, I've seen it zero times
	alt_num = 0
	#before I've seen anything, I don't have an alternative allele
	alt_allele = ""
	#loop over all the possible alternative alleles
	for allele in alt_alleles:
		if alt_alleles[allele] > alt_num:
			#if I've seen this allele more than any other, this is my alternative allele
			alt_num = alt_alleles[allele]
			alt_allele = allele
	if alt_num > 0:
		#if there are more 0 alternative reads, then this is a potentially variable site, so compute GLs and print
		GLref, GLalt = compute_genotype_likelihood(ref_num,alt_num)
		outFile.write("%s\t%d\t%s\t%s\t%f\t%f\n"%(chrom,pos,ref_allele,alt_allele,GLref,GLalt))
	
