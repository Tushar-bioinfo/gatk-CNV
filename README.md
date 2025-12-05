# GATK4-CNV-Pipeline  
Copy Number Variation Analysis for Cancer RNA-Seq Data

This project implements a reproducible workflow for computing copy number variations (CNVs) from large-scale dbGaP/GDC RNA-seq cohorts using GATK4. The pipeline is organized into three major stages:  
(1) **Pre-processing** (SRA → FASTQ → BAM),  
(2) **CNV Processing** (GATK CNV tools), and  
(3) **Post-processing** (gene-level CNV matrix + survival analysis).

All steps are optimized for HPC execution on Slurm systems and support downstream statistical or clinical analyses.

---

## Project Structure  

```
├── workflow/
│   ├── 01-pre_processing/          # SRA → FASTQ → BAM preparation
│   │   ├── 01_pre.sh
│   │   ├── 02_pre.sh
│   │   ├── 03_pre.py
│   │   └── 04_pre.ipynb
│   │
│   ├── 02-processing/              # GATK4 CNV workflow
│   │   ├── 01_main.sh
│   │   ├── 02_main.sh
│   │   ├── 03_main.sh
│   │   ├── 04_main.sh
│   │   ├── 05_main.sh
│   │   ├── 06_main.sh
│   │   ├── 07_main.sh
│   │   └── 08_main.sh
│   │
│   └── 03-post_processing/         # Gene-level CNV + survival analysis
│       ├── main_01.R
│       ├── main_02.R
│       └── survival_analysis.ipynb
│
└── README.md
```
---

## Pre-Processing

The pre-processing stage converts raw SRA files into aligned, analysis-ready BAM files.

Steps include:

- Download and verify SRA accessions  
- Extract FASTQ files using `fasterq-dump` (scratch-based execution)  
- Align FASTQ with BWA-mem2  
- Sort and index BAM files using samtools  
- Perform QC and metadata preparation (Python + Jupyter)

Location: `workflow/01-pre_processing/`

---

## CNV Processing (GATK4)

This stage runs the full GATK CNV workflow:

- CollectReadCounts for tumor and normal samples  
- Build Panel of Normals (CreateReadCountPanelOfNormals)  
- Apply DenoiseReadCounts using the PoN  
- Segment genomes using ModelSegments  
- Export CNV segments and denoised copy-ratio files  

Location: `workflow/02-processing/`

---

## Post-Processing & Analysis

After segmentation, CNV results are aggregated and prepared for downstream analyses.

Includes:

- Mapping CNV segments to gene coordinates (R)  
- Creating a unified gene-level CNV matrix  
- Performing survival analysis (Kaplan–Meier and Cox models) in Jupyter notebooks  

Location: `workflow/03-post_processing/`

---

## Requirements

### Software
- GATK 4.x  
- BWA-mem2  
- samtools  
- SRA Toolkit  
- Python 3.9+  
- R 4.0+  
- Slurm HPC scheduler  

### R Packages
- dplyr  
- data.table  
- GenomicRanges  
- readr  

### Python Packages
- pandas  
- numpy  

---

### Outputs
- Denoised copy-ratio files
- CNV segmentation tables
- Gene-level CNV matrix
- Survival analysis plots and statistical summaries
