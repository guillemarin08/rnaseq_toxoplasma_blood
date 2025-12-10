# RNA-seq Analysis — *Toxoplasma gondii* Infection in Mouse Blood

**Author:** Guillermo Marín García  
**Email:** guillermo.maringarcia@students.unibe.ch  
**Course:** RNA-seq Analysis, University of Bern  
**Dataset:** GEO GSE119855 (Singhania et al., 2019)

---

## Project Overview
This repository contains the scripts used to process and analyze RNA-seq data from blood samples of *Mus musculus* under *Toxoplasma gondii* infection.  

The dataset includes:
- **Genotypes:** Wild-type (WT) and interferon receptor double knockout (DKO; *Ifnar1−/− Ifngr1−/−*)
- **Conditions:** Infected (*Case*) and uninfected (*Control*)
- **Sequencing:** Paired-end Illumina HiSeq 4000

The analysis follows a standard bulk RNA-seq workflow from raw FASTQ files to differential expression and GO enrichment.

---

## Sample Information

| Sample ID     | Genotype | Condition |
|---------------|----------|-----------|
| SRR7821949    | WT       | Case      |
| SRR7821950    | WT       | Case      |
| SRR7821951    | WT       | Case      |
| SRR7821952    | WT       | Case      |
| SRR7821953    | WT       | Case      |
| SRR7821968    | WT       | Control   |
| SRR7821969    | WT       | Control   |
| SRR7821970    | WT       | Control   |
| SRR7821954    | DKO      | Case      |
| SRR7821955    | DKO      | Case      |
| SRR7821956    | DKO      | Case      |
| SRR7821957    | DKO      | Case      |
| SRR7821971    | DKO      | Control   |
| SRR7821972    | DKO      | Control   |
| SRR7821973    | DKO      | Control   |

---

## Workflow Summary
1. **Quality control of raw reads** — FastQC  
2. **Trimming and filtering** — fastp  
3. **Post-trim QC aggregation** — MultiQC  
4. **Genome indexing and alignment** — HISAT2  
5. **BAM processing** — Samtools  
6. **Read quantification** — featureCounts  
7. **Differential expression analysis** — DESeq2  
8. **GO enrichment** — clusterProfiler  

All steps are implemented through Bash scripts and an R script for downstream analysis.

---

## Repository Structure

.
├── scripts/ # All analysis scripts (Bash + R)
│ ├── 01_fastqc_raw.sh
│ ├── 02_fastp.sh
│ ├── 03_fastqc_fastp.sh
│ ├── 04_multiqc_fastp.sh
│ ├── 05_multiqc_raw_vs_fastp.sh
│ ├── 06_hisat2_index.sh
│ ├── 07_hisat2_mapping.sh
│ ├── 08_featurecounts.sh
│ ├── 09_final_multiqc.sh
│ └── 10_DESeq2_analysis.R
├── README.md # Project documentation
└── .gitignore # Excludes large/unnecessary files

---

## How to Run the Workflow

Each step can be submitted to the HPC cluster via Slurm:

sbatch scripts/01_fastqc_raw.sh
sbatch scripts/02_fastp.sh
sbatch scripts/06_hisat2_index.sh
sbatch scripts/07_hisat2_mapping.sh
sbatch scripts/08_featurecounts.sh


## DESeq2 and enrichment analysis are run in R:

Rscript scripts/10_DESeq2_analysis.R



## Software Versions
FastQC v0.12.1
fastp v0.24.1
MultiQC v1.19
HISAT2 v2.2.1
Samtools v1.20
Subread/featureCounts v2.0.6
R v4.3.2
DESeq2 v1.42.1
clusterProfiler v4.10.1