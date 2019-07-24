# This script parses the VCF file.  It counts the total number of mutations (with less than 5% prevalence)
#!/usr/bin/env python 

from argparse import ArgumentParser
vcfFile = False
bamReadcount = False
import math


parser = ArgumentParser()
parser.add_argument("-v", "--vcf", dest='vcfFile',
                    help="vcf file", metavar="FILE")
parser.add_argument("-r", "--bam_readcount", dest='bamReadcount',
                    help="reference genome in FASTA format", metavar="FILE")
parser.add_argument("-o", "--output", dest='output',
                    help="output filename", metavar="STR")

args = parser.parse_args()
vcfFile = args.vcfFile
bamReadcount = args.bamReadcount
output = args.output


#First, get the number of mutations by parsing the VCF file.

f = open(vcfFile,"r")
g = open(output,"w")
mutationCount = 0
for line in f:
    if "#" in line:
    	pass
    else:
    	t = line.split()
    	info = t[7].split(";") 
    	frequency = info[5].split("=")
    	frequency = frequency[1].split(",")
    	mutFreq = 0
    	for i in frequency:
    		mutFreq += int(i)
    	depth = info[7].split("=")
    	depth = int(depth[1])
    	if mutFreq / depth < 0.05:
    		mutationCount += mutFreq

		
f.close()


#Next, parse bam-readcount to find the total number of bases

f = open(bamReadcount,"r")
baseCount = 0
for line in f:
	t = line.split()
	baseCount += int(t[3])
f.close()

#write the total mutations and the total bases sequenced to file, and calculate the percentage of bases that have mutations

g.write("total mutations = " + str(mutationCount) + "\n")
g.write("total bases = " + str(baseCount)+ "\n")
incorrect = mutationCount/baseCount
g.write("Mutation rate = " + str(incorrect*100)+"%\n")
#calculate what phred would be for this error rate
if incorrect !=0:
	phred = math.log(incorrect,10) * -10
	g.write("Phred = " + str(phred))
else:
	g.write("No mutations, phred cannot be calculated")
g.close()
