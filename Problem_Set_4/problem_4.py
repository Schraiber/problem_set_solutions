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

def compute_genotype_likelihood(n, alt_num,e=0.01):
	ref_num = n - alt_num
	GLref = binom(n,alt_num)*e**alt_num*(1.0-e)**ref_num
	GLalt = binom(n,alt_num)*(1.0-e)**alt_num*e**ref_num
	return GLref, GLalt

bamFile = pysam.AlignmentFile("BC187_mapped.sorted.bam")

fastaFile = open("yeast_genome.fa")

refGenome = read_fasta(fastaFile)

outFile = open("yeast_genotype_likelihoods.txt","w")

for pileupcolumn in bamFile.pileup():
	chrom = pileupcolumn.reference_name
	pos = pileupcolumn.pos
	n = pileupcolumn.n
	ref_allele = refGenome[chrom][pos]
	alt_alleles = {}
	for pileupread in pileupcolumn.pileups:
		if pileupread.is_del or pileupread.is_refskip: continue
		read = pileupread.alignment.query_sequence[pileupread.query_position]
		if read != ref_allele:
			if read in alt_alleles:
				alt_alleles[read] += 1
			else:
				alt_alleles[read] = 0
	alt_num = 0
	alt_allele = ""
	for allele in alt_alleles:
		if alt_alleles[allele] > alt_num:
			alt_num = alt_alleles[allele]
			alt_allele = allele
	if alt_num > 0:
		GLref, GLalt = compute_genotype_likelihood(n,alt_num)
		outFile.write("%s\t%d\t%s\t%s\t%f\t%f\n"%(chrom,pos,ref_allele,alt_allele,GLref,GLalt))
	
