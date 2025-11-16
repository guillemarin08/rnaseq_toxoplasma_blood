# RNA-seq Analysis: Toxoplasma Infection in Mouse Blood

Author: **Guillermo Marín García**  
Email: **guillermo.maringarcia@students.unibe.ch**  
GitHub: **https://github.com/guillemarin08**

This project is part of the RNA-seq course at the University of Bern.  
It analyzes Toxoplasma gondii infection using mouse blood RNA-seq data 
from Singhania et al. (2019).

Tissue: Blood
Table describing experimental groups for each sample. Samples are either from wildtype (WT) mice or from double-knockout (DKO) mice. 
Since the reads were produced in paired-end mode, there are 2 files per sample, with read 1 and read 2 respectively.

Sample = ID as it is in the fastq file name
Case = Infected
Control = Uninfected
WT = Wildtype
DKO = Double-Knockout

Sample	Group
SRR7821949	Blood_WT_Case
SRR7821950	Blood_WT_Case
SRR7821951	Blood_WT_Case
SRR7821952	Blood_WT_Case
SRR7821953	Blood_WT_Case
SRR7821968	Blood_WT_Control
SRR7821969	Blood_WT_Control
SRR7821970	Blood_WT_Control
SRR7821954	Blood_DKO_Case
SRR7821955	Blood_DKO_Case
SRR7821956	Blood_DKO_Case
SRR7821957	Blood_DKO_Case
SRR7821971	Blood_DKO_Control
SRR7821972	Blood_DKO_Control
SRR7821973	Blood_DKO_Control

## Workflow
1. Raw data QC — FastQC
2. Read trimming — fastp/Trimmomatic
3. Post-trim QC — FastQC + MultiQC
4. Alignment — HISAT2
5. Sorting & indexing — samtools
6. Counting — featureCounts
7. DESeq2 analysis
8. GO enrichment


## Structure

