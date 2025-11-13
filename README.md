# RNA-seq Analysis: Toxoplasma Infection in Mouse Blood

Author: **Guillermo Marín García**  
Email: **guillermo.maringarcia@students.unibe.ch**  
GitHub: **https://github.com/guillemarin08**

This project is part of the RNA-seq course at the University of Bern.  
It analyzes Toxoplasma gondii infection using mouse blood RNA-seq data 
from Singhania et al. (2019).

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

