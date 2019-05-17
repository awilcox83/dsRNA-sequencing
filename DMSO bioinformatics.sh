# Shell script to run analyses associated with the paper. See README file.  
# Download all FASTQ files to the same directory. The following bash scripts will carry out our 
# analyses on every pair of FASTQ files in this folder.

#Required Packages and version tested
#sickle 1.33
#bowtie2 2.2.6
#pear 0.9.6
#samtools 1.9
#bam-readcountunstable (commit nogit)

#Example of installation on required packages in each of 
#Install required packages
#module load sickle/1.33
#module load bowtie2/2.2.6
#conda install -c bioconda pear #Tested in version v0.9.6
#conda install -c bioconda samtools #Tested in version 1.9

#Check on each program's installation
cutadapt -h
sickle -h
pear
samtools --help 
bowtie2 -h
bam-readcount -v

# Download and rename the sequence files from SRA via EBI

#DMSO Concentration experiments with influenza virus A/California/07/2009 and bacteriophage phi-6 (experiment 1)
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR872/009/SRR8729329/SRR8729329_1.fastq.gz -O influenza_0_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR872/000/SRR8729330/SRR8729330_1.fastq.gz -O influenza_15_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR872/001/SRR8729331/SRR8729331_1.fastq.gz -O influenza_50_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR872/002/SRR8729332/SRR8729332_1.fastq.gz -O influenza_90_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR872/009/SRR8729329/SRR8729329_2.fastq.gz -O influenza_0_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR872/000/SRR8729330/SRR8729330_2.fastq.gz -O influenza_15_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR872/001/SRR8729331/SRR8729331_2.fastq.gz -O influenza_50_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR872/002/SRR8729332/SRR8729332_2.fastq.gz -O influenza_90_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR872/000/SRR8729120/SRR8729120_1.fastq.gz -O phi6_0_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR872/009/SRR8729119/SRR8729119_1.fastq.gz -O phi6_15_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR872/008/SRR8729118/SRR8729118_1.fastq.gz -O phi6_50_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR872/007/SRR8729117/SRR8729117_1.fastq.gz -O phi6_90_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR872/000/SRR8729120/SRR8729120_2.fastq.gz -O phi6_0_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR872/009/SRR8729119/SRR8729119_2.fastq.gz -O phi6_15_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR872/008/SRR8729118/SRR8729118_2.fastq.gz -O phi6_50_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR872/007/SRR8729117/SRR8729117_2.fastq.gz -O phi6_90_R2.fastq.gz

#Experiments to test the effect of nuclease treatment, lysate concentration, and DMSO on viral read recovery (experiment 2)
# In the following filenames, N = treated with nucleases, C = concentrated sample, D = treated with DMSO
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR903/001/SRR9030131/SRR9030131_1.fastq.gz -O Phi6_NCD_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR903/001/SRR9030131/SRR9030131_2.fastq.gz -O Phi6_NCD_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR903/002/SRR9030132/SRR9030132_1.fastq.gz -O Phi6_NC_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR903/002/SRR9030132/SRR9030132_2.fastq.gz -O Phi6_NC_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR903/003/SRR9030133/SRR9030133_1.fastq.gz -O Phi6_N_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR903/003/SRR9030133/SRR9030133_2.fastq.gz -O Phi6_N_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR903/004/SRR9030134/SRR9030134_1.fastq.gz -O Phi6_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR903/004/SRR9030134/SRR9030134_2.fastq.gz -O Phi6_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR903/005/SRR9030135/SRR9030135_1.fastq.gz -O Phi6_C_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR903/005/SRR9030135/SRR9030135_2.fastq.gz -O Phi6_C_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR903/006/SRR9030136/SRR9030136_1.fastq.gz -O Phi6_ND_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR903/006/SRR9030136/SRR9030136_2.fastq.gz -O Phi6_ND_R2.fastq.gz

#These files were analyzed with a custom virus discovery pipeline, so commenting out from further analysis. See paper for details
# wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR872/005/SRR8729125/SRR8729125_1.fastq.gz -O cocktail_0_R1.fastq.gz
# wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR872/004/SRR8729124/SRR8729124_1.fastq.gz -O cocktail_15_R1.fastq.gz
# wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR872/003/SRR8729123/SRR8729123_1.fastq.gz -O cocktail_50_R1.fastq.gz
# wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR872/002/SRR8729122/SRR8729122_1.fastq.gz -O cocktail_90_R1.fastq.gz
# wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR872/005/SRR8729125/SRR8729125_2.fastq.gz -O cocktail_0_R2.fastq.gz
# wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR872/004/SRR8729124/SRR8729124_2.fastq.gz -O cocktail_15_R2.fastq.gz
# wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR872/003/SRR8729123/SRR8729123_2.fastq.gz -O cocktail_50_R2.fastq.gz
# wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR872/002/SRR8729122/SRR8729122_2.fastq.gz -O cocktail_90_R2.fastq.gz

#Number of reads per file
wc -l *.fastq.gz | awk '{print $1/4 " " $2}'

# cutadapt was used to trim adapters and primers from reads. Primer and adapter sequences for Illumina Nextera library prep kit retrieved from https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences-1000000002694-10.pdf“.
for file in *_R1.fastq.gz
do
	prefix=${file%_R1.fastq.gz}
	cutadapt ${prefix}_R1.fastq.gz ${prefix}_R2.fastq.gz -m 25 -n 2 -o ${prefix}_R1.trimmed.fastq.gz -p ${prefix}_R2.trimmed.fastq.gz \
	-b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
	-a CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGG -a CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGG \
	-a CAAGCAGAAGACGGCATACGAGATCTAGTACGGTCTCGTGGGCTCGG -a CAAGCAGAAGACGGCATACGAGATTTCTGCCTGTCTCGTGGGCTCGG \
	-a CAAGCAGAAGACGGCATACGAGATGCTCAGGAGTCTCGTGGGCTCGG -a CAAGCAGAAGACGGCATACGAGATAGGAGTCCGTCTCGTGGGCTCGG \
	-a CAAGCAGAAGACGGCATACGAGATCATGCCTAGTCTCGTGGGCTCGG -a AATGATACGGCGACCACCGAGATCTACACTAGATCGCTCGTCGGCAGCGTC \
	-a AATGATACGGCGACCACCGAGATCTACACCTCTCTATTCGTCGGCAGCGTC > ${prefix}.txt
done
#Results of cutadapt command are in .txt files for each file (same basename, i.e. influenza_0.txt)

# Sickle was used to trim for quality
for file in *_R1.trimmed.fastq.gz
do
	prefix=${file%_R1.trimmed.fastq.gz}
	sickle pe -f ${prefix}_R1.trimmed.fastq.gz -r ${prefix}_R2.trimmed.fastq.gz -t sanger -o ${prefix}_R1.qc.fastq.gz -p ${prefix}_R2.qc.fastq.gz -s ${prefix}.qcsingles.fastq.gz -n > sickle_${prefix}.txt
done
#Results of sickle command are in .txt files for each file (same basename, i.e. influenza_0.txt)

# PEAR was used to merge paired end reads.  This command was repeated for each set of paired fastq files
for file in *_R1.qc.fastq.gz
do
	prefix=${file%_R1.qc.fastq.gz}
	PEAR -f ${prefix}_R1.qc.fastq.gz -r ${prefix}_R2.qc.fastq.gz -n 25 -o ${prefix}.merged > pear_${prefix}.txt
done
#Results of pear command are in .txt files for each file (same basename, i.e. influenza_0.txt)

# bowtie2 was used to index the reference genomes
bowtie2-build influenza.fasta influenza_ref
bowtie2-build phi6.fasta phi6_ref

#Check in on number of reads after all QC Steps
#Number of reads in unassembled and merged/assembled FASTQ's per file
wc -l *.fastq | awk '{print $1/4 " " $2}'

#Number of reads only in merged/assembled FASTQ's, used going forward
wc -l *.merged.assembled.fastq | awk '{print $1/4 " " $2}'

# reads were mapped to the reference genomes with bowtie 2.  Note that this only maps reads from experiment 1.  Experiment 2 reads will be mapped separately.
for file in influenza*.merged.assembled.fastq
do
	prefix=${file%.merged.assembled.fastq}
	bowtie2 -x influenza_ref -U ${prefix}.merged.assembled.fastq -S ${prefix}.mapped.SAM 2> bowtie_${prefix}.txt 
done
#Results of bowtie2 command are in .txt files for each file (same basename, i.e. influenza_0.txt)

for file in phi6*.merged.assembled.fastq
do
	prefix=${file%.merged.assembled.fastq}
	bowtie2 -x phi6_ref -U ${prefix}.merged.assembled.fastq -S ${prefix}.mapped.SAM  2> bowtie_${prefix}.txt
done
#Results of bowtie2 command are in .txt files for each file (same basename, i.e. influenza_0.txt)

# SAM files were converted to BAM, sorted, and indexed
for file in *.mapped.SAM
do
	prefix=${file%.mapped.SAM}
	samtools view -bS ${prefix}.mapped.SAM > ${prefix}.mapped.BAM
	samtools sort ${prefix}.mapped.BAM -o ${prefix}.mapped.sorted.BAM
	samtools index ${prefix}.mapped.sorted.BAM
done

# bam-readcount was used to determine the read depth at each position
# NOTE: Known issue with warnings, can disregard
for file in *.mapped.sorted.BAM
do
	prefix=${file%.mapped.sorted.BAM}
	bam-readcount -b 20 -w 1 ${prefix}.mapped.sorted.BAM > ${prefix}.tab
done

#This final readcount file tab-delimited (*.tab) is used going forward to generate plots in R


# reads from experiment 2 mapped to phi6 reference genome.  
for file in Phi6*.merged.assembled.fastq
do
	prefix=${file%.merged.assembled.fastq}
	bowtie2 -x phi6_ref -U ${prefix}.merged.assembled.fastq -S ${prefix}.mapped.SAM  2> bowtie_${prefix}.txt
done

# Mapping data from bowtie converted to tab-delimited file file for use in R
printf 'Name\tTotal\tUnpaired\tUnaligned\tAlignedOnce\tAlignedMultiple\tPercentAligned\n' >> mapping_percentages.tsv
for file in bowtie_Phi6*
do
	prefix=${file%.txt}
	title=${prefix#"bowtie_"}
	printf ${title} >> mapping_percentages.tsv
	printf '\t' >> mapping_percentages.tsv
	while read word _; do printf '%s\t' "$word"; done < ${prefix}.txt >> mapping_percentages.tsv
	printf '\n' >> 	mapping_percentages.tsv
done


#Run R scripts to generate coverage plots and CSV files of data
Rscript read_coverage_plots.R
Rscript mapping_percentage_plots.R
