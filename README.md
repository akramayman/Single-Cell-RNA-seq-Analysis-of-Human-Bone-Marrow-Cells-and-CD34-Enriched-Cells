# 🧬 scRNA-seq Analysis of Human Bone Marrow and CD34+ Enriched Cells
![R Package](https://img.shields.io/badge/R%20Package-Seurat-blue)
![R Package](https://img.shields.io/badge/R%20Package-CellChat-green)
![R Package](https://img.shields.io/badge/R%20Package-Monocle3-red)
![R Package](https://img.shields.io/badge/R%20Package-ggplot2-yellow)
 

This study analyzes scRNA-seq data from human bone marrow mononuclear cells (BMMC) and CD34+ hematopoietic progenitor cells collected from multiple donors. The goal is to reconstruct the cellular landscape of the human bone marrow and identify major immune and stem cell populations through transcriptomic profiling.

The analysis follows a complete standard single-cell workflow, implemented primarily in R using Seurat, including quality control, batch correction, dimensionality reduction, clustering, and cell type annotation.

---

## 🧪  Table of Contents

- [Dataset Overview](#dataset-overview)
- [Quality Control](#quality-control)
- [Batch Correction](#batch-correction)
- [Dimensionality Reduction](#dimensionality-reduction)
- [Clustering](#clustering)
- [Cell Type Annotation](#cell-type-annotation)
- [Differential Expression Analysis](#differential-expression-analysis)

---

## 📦 Dataset Overview

The dataset consists of four samples from two cell populations across multiple donors:

| Sample      | Cell Type | Number of Cells | Number of Genes |
|-------------|-----------|-----------------|-----------------|
| BMMC-D1T1   | BMMC      | 6,270           | 14,065          |
| BMMC-D1T2   | BMMC      | 6,332           | 14,084          |
| CD34-D2T1   | CD34+     | 2,424           | 13,441          |
| CD34-D3T1   | CD34+     | 5,752           | 13,432          |

- **BMMC**: Bone Marrow Mononuclear Cells (mature immune cells)
- **CD34+**: Hematopoietic stem and progenitor cells enriched for CD34 surface marker

The dataset used in this project can be downloaded here:

- [Download scRNA-seq bone marrow dataset (ZIP)](https://icbb-share.s3.eu-central-1.amazonaws.com/single-cell-bioinformatics/scbi_ds1.zip)
  
---

## 🔬 Quality Control

Violin plots (VlnPlot) were generated for each sample to assess three key QC metrics:

- **nFeature_RNA** — number of unique genes detected per cell
- **nCount_RNA** — total UMI counts per cell
- **percent.mt** — percentage of mitochondrial gene expression

Key observations:
- The BMMC samples (D1T1, D1T2) show a sharp, narrow distribution in both nFeature and nCount, consistent with a heterogeneous but well-captured mature cell population.
- The CD34 samples show a broader distribution, reflecting the expected transcriptional diversity of progenitor cells.
- Mitochondrial content (percent.mt) is essentially zero across all samples, indicating high cell viability and minimal apoptotic contamination.

---

## ⚖️ Batch Correction

UMAP plots were generated before and after batch correction, colored by several metadata variables to identify sources of technical variation:

| Variable     | Observation |
|--------------|-------------|
| `orig.ident` | BMMC and CD34 cells were initially separated; after correction they integrate while preserving biological identity |
| `Donor`      | D1 (BMMC) dominates numerically; D2 and D3 (CD34) mix well post-correction |
| `Sex`        | Predominantly female (F) donors; male (M) cells distributed across shared clusters after correction |
| `Replicate`  | T1/T2 replicates mix well after correction, confirming technical reproducibility |

Batch correction was successful in removing technical donor/replicate effects while preserving biologically meaningful cell type structure.

---

## 📉 Dimensionality Reduction

### PCA

A scree plot was used to select the number of principal components (PCs) for downstream analysis. The "elbow" occurs around **PC 6–7**, where explained variance begins to plateau. The first **10 PCs** were selected as they capture the majority of the biological signal.

### UMAP

UMAP was run on the top 10 PCs, revealing **15 distinct clusters** (0–14) with clear separation between major hematopoietic lineages. The 2D embedding shows:
- A large central myeloid/progenitor compartment
- Distinct lymphoid populations
- Isolated erythroid and NK/T cell clusters

---

## 🧩 Clustering

Seurat graph-based clustering identified **15 clusters** (0–14) at the chosen resolution. Side-by-side UMAP plots colored by `orig.ident` and `seurat_clusters` confirm that:
- CD34+ progenitor cells predominantly occupy the upper clusters (enriched in stem/progenitor identities)
- BMMC cells span the full landscape, including mature myeloid, lymphoid, and erythroid compartments

---

## Cell Type Annotation

### Automatic Annotation

Automated annotation assigned the following cell type labels:

`B_cell`, `BM`, `BM & Prog.`, `CMP`, `DC`, `Erythroblast`, `GMP`, `HSC_-G-CSF`, `MEP`, `Monocyte`, `Myelocyte`, `Neutrophils`, `NK_cell`, `Platelets`, `Pre-B_cell_CD34-`, `Pro-B_cell_CD34+`, `Pro-Myelocyte`, `T_cells`, `Tissue_stem_cells`, `NA`

### Manual Annotation

Manual annotation using canonical marker genes produced the following labels:

`T-Cells`, `CD14` (Monocytes), `GMP`, `Plasma`, `Basophils`, `cDC`, `Pre-Cell`, `Erythrocytes`, `CLP`, `CD8`, `B-Cell`, `CD4`, `LMPP`, `NK`, `CD16`

### Comparison

Both approaches yielded broadly similar cluster distributions. The key difference is granularity: automatic annotation treated the T cell population as a single cluster, while manual annotation resolved it into three distinct subsets — **T-Cells**, **CD4**, and **CD8** — providing finer immunological resolution. This means the automatic approach annotates more cells per label, while manual annotation offers greater biological specificity.

---

## 📊 Differential Expression Analysis

Marker gene expression was visualized using violin plots and UMAP feature plots for three canonical lineage markers:

| Gene   | Expression Pattern |
|--------|--------------------|
| **CD19**  | Highly enriched in clusters 6, 10, and 12 — consistent with B cell identity |
| **CD3D**  | Strongly expressed in clusters 0, 9, and 11 — confirming T cell populations |
| **NKG7**  | Highest in clusters 11 and 13 — marking NK and cytotoxic T cells |

### Pairwise Comparisons (Volcano Plots)

**B-Cells vs T-Cells**
- Large number of significantly differentially expressed genes in both directions
- Genes upregulated in B cells include B cell receptor components; genes upregulated in T cells include T cell receptor and co-stimulatory molecules

**CD4 T Cells vs CD14 Monocytes**
- Strong transcriptional separation between the two lineages
- Monocyte-enriched genes (upregulated in CD14) reflect innate immune and phagocytic functions; T cell genes reflect adaptive immune signaling

---

## Tools & Technologies

- **R / Seurat** — QC, normalization, clustering, annotation
- **UMAP** — non-linear dimensionality reduction
- **PCA** — linear dimensionality reduction and PC selection
- **Harmony / integration** — batch correction
- **ggplot2** — visualization

--- 

## Citation 

If you use this repository, please cite:
Akram Abushmais.
