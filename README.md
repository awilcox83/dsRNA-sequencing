# dsRNAseq: Code for "Next-generation sequencing of double stranded RNA is greatly improved by treatment with the inexpensive denaturing reagent DMSO"

This github repository includes code (and links to data) from the manuscript: 
Next-generation sequencing of double stranded RNA is greatly improved by treatment with the inexpensive denaturing reagent DMSO
Alexander H. Wilcox, Eric Delwart, and Samuel L. Díaz Muñoz

If you are reading or using this, let us know how these data were useful for you. Always open to collaborate! Please contact us!

If you use these data and code, please cite the bioRxiv preprint or paper.

## Quick Start
1. Make sure packages are installed (see comments in DMSO_bioinformatics.sh)
2. chmod +x DMSO_bioinformatics.sh
3. ./DMSO_bioinformatics.sh

### 1. Project Description

The origin of this study was in our efforts to sequence double-stranded RNA phages in the Cystoviridae (phi6 and its relatives) in a cost-effective and high-throughput manner. Alex hit upon earlier literature using DMSO to improve reverse transcriptase PCR (RT_PCR) and Sanger Sequencing, and tried using DMSO to imporve next-generation sequencing. This manuscript details our effort to validate this methodology. We found that DMSO greatly improved sequencing read recovery even with very low starting nucleic acid concentrations. We also tested whether DMSO affected read recovery from ssRNA viral genomes (Influenza), and whether it improved recovery from metagenomic samples - see abstract below for details. This repository has the code underlying the results presented in the manuscript.

Abstract:
Double stranded RNA (dsRNA) is the genetic material of important viruses and a key component of RNA interference-based immunity in eukaryotes. Previous studies have noted difficulties in determining the sequence of dsRNA molecules that have affected studies of immune function and estimates of viral diversity in nature. Dimethyl sulfoxide (DMSO) has been used to denature dsRNA prior to the reverse transcription stage to improve RT-PCR and Sanger sequencing. We systematically tested the utility of DMSO to improve sequencing yield of a dsRNA virus (Φ6) in a short-read next generation sequencing platform. DMSO treatment improved sequencing read recovery by over two orders of magnitude, even when RNA and cDNA concentrations were below the limit of detection. We also tested the effects of DMSO on a mock eukaryotic viral community and found that dsRNA virus reads increased with DMSO treatment. Furthermore, we provide evidence that DMSO treatment does not adversely affect recovery of reads from a single-stranded RNA viral genome (Influenza A/California/07/2009). We suggest that up to 50% DMSO treatment be used prior to cDNA synthesis when samples of interest are composed of or may contain dsRNA.

### 2. Data and Reference Files
Raw data consists of sequencing output from the Illumina MiSeq platform, i.e. FASTQ files. These files are available in the Sequencing Read Archive (SRA, NCBI) under three BioProject Accessions:  

Bacteriophage phi-6 (Cystoviridae):
Accession: PRJNA527098 (ID: 527098) https://www.ncbi.nlm.nih.gov/bioproject/527098

Influenza Virus (A/California/07/2009):
Accession: PRJNA527101 (ID: 527101) https://www.ncbi.nlm.nih.gov/bioproject/PRJNA527101

Mock Eukaryotic Viral Community (NIBSC reagent 11/242-001)
Accession: PRJNA527100 (ID: 527100)
https://www.ncbi.nlm.nih.gov/bioproject/527100 

We also include the FASTA files that serve as reference sequences for influenza (influenza.fasta) and phi6 (phi6.fasta).

### 3. Code
Below are descriptions of the code files used to generate the tables, figures, and statistics in the paper.

1) DMSO_bioinformatics.sh: This file is shell script that downloads and processes raw sequencing reads from the BioProjects above. The commands process and map sequence data to reference genomes and culminate in files determining read depth at each position in the phi-6 and influenza virus genomes. Please note the comments at the beginning of the file, which outline the software requirements.  This script calls the following two R scripts.

2) read_coverage_plots.R: This file is the R script used to generate the figures showing the coverage per base. It writes data used as CSV, and saves plots as PDF.

3) mapping_percentages.R: This file is the R script used to generate the figure showing the percentage of reads that map to phi-6 under different treatments before sequencing.  It saves its plot as a PDF.
