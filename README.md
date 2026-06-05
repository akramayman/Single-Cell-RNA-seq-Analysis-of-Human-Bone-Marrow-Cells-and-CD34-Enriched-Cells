# Single-Cell-RNA-seq-Analysis-of-Human-Bone-Marrow-Cells-and-CD34-Enriched-Cells

Single-cell RNA sequencing (scRNA-seq) has become a powerful tool for studying cellular heterogeneity and identifying distinct cell populations within complex biological systems. This project presents a comprehensive scRNA-seq analysis of human bone marrow mononuclear cells (BMMCs) and CD34+ hematopoietic progenitor cells obtained from multiple donors and biological replicates.

The primary objective of this study is to characterize the cellular composition of hematopoietic populations by applying a standard single-cell analysis workflow, including quality control, normalization, dimensionality reduction, clustering, cell-type annotation, and differential gene expression analysis. To ensure data quality and reliability, doublets were identified and removed using DoubletFinder, while batch effects arising from donor and replicate differences were corrected through Seurat's integration framework.

Both automated annotation using SingleR and manual annotation based on known marker genes were performed to identify major immune and progenitor cell populations. Additionally, differential expression analyses were conducted to investigate transcriptional differences between selected cell types, providing insights into lineage-specific gene expression patterns.

This project demonstrates an end-to-end scRNA-seq analysis pipeline using the Seurat ecosystem in R and serves as a practical example of single-cell bioinformatics workflows for studying hematopoietic cell diversity and gene expression dynamics.


## Dataset

### Source

Dataset download link:

https://icbb-share.s3.eu-central-1.amazonaws.com/single-cell-bioinformatics/scbi_ds1.zip

### Samples

| Sample | Description |
|----------|----------|
| BMMC-D1T1 | Bone Marrow Mononuclear Cells - Donor 1 Replicate 1 |
| BMMC-D1T2 | Bone Marrow Mononuclear Cells - Donor 1 Replicate 2 |
| CD34-D2T1 | CD34+ Cells - Donor 2 Replicate 1 |
| CD34-D3T1 | CD34+ Cells - Donor 3 Replicate 1 |

### Dataset Statistics

| Sample | Cells | Genes |
|----------|-------|-------|
| BMMC-D1T1 | 6270 | 14065 |
| BMMC-D1T2 | 6332 | 14084 |
| CD34-D2T1 | 2424 | 13441 |
| CD34-D3T1 | 5752 | 13432 |

---

## Installation

### Clone Repository

git clone https://github.com/username/repository.git

### Create Environment

conda env create -f environment.yaml
conda activate scrnaseq

---

## Analysis Workflow

1. Data Loading
2. Quality Control
3. Normalization
4. Feature Selection
5. PCA
6. Doublet Detection (DoubletFinder)
7. Data Integration & Batch Correction
8. Clustering
9. Automatic Cell Annotation (SingleR)
10. Manual Cell Annotation
11. Differential Expression Analysis

---

## Main Packages

- Seurat
- DoubletFinder
- SingleR
- celldex
- Monocle3
- CellChat
- tidyverse
- ggpubr
