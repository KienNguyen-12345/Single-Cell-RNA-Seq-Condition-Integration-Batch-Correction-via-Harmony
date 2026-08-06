# Single-Cell RNA-Seq Condition Integration & Batch Correction via Harmony

## Overview
This repository contains an R-based single-cell RNA sequencing (scRNA-seq) workflow demonstrating dataset integration across experimental conditions using the **Harmony** algorithm within the **Seurat** framework. The analysis integrates peripheral blood mononuclear cells (PBMCs) from control and interferon-beta (IFN-β) stimulated conditions to correct condition-level batch variation and enable cell-type-specific differential evaluation.

---

## Dataset
* **Source:** `ifnb` dataset from `SeuratData`
* **Study Design:** PBMCs from 8 lupus patients split into two experimental conditions:
  * **CTRL:** Control / Untreated
  * **STIM:** Simulated / IFN-β treated
* **Goal:** Overlay shared immune cell populations across conditions to enable unbiased comparative transcriptomic analysis.

---

## Analytical Pipeline & Workflow

### 1. Quality Control & Filtering
* Evaluated mitochondrial read percentages (`pattern = "^MT-"`) and UMI counts per cell.
* Filtered low-quality cells maintaining `nCount_RNA > 800`.

### 2. Standard Pre-processing
* Log-normalization (`NormalizeData`) and highly variable feature selection (`FindVariableFeatures`).
* Data scaling (`ScaleData`) and Principal Component Analysis (`RunPCA`).
* Baseline UMAP reduction showing clear condition-driven separation prior to batch correction.

### 3. Automated Cell-Type Annotation
* Converted Seurat object to `SingleCellExperiment` format.
* Mapped cell identities using `SingleR` against the `celldex::HumanPrimaryCellAtlasData` reference dataset.
* Added predicted cell labels back to Seurat metadata for downstream verification.

### 4. Integration with Harmony
* Executed `RunHarmony(group.by.vars = 'stim')` on low-dimensional PCA embeddings.
* Performed downstream graph-based clustering (`FindNeighbors`, `FindClusters`) and UMAP reduction using the corrected `harmony` embeddings rather than raw PCs.

### 5. Visualization & Validation
* Generated side-by-side comparative UMAP plots (`before` vs. `after` Harmony integration).
* Verified successful alignment and overlay of corresponding immune cell lineages across `CTRL` and `STIM` groups.

---

## Key Dependencies & R Packages
* **Seurat / SeuratData:** Core scRNA-seq processing & data container
* **harmony:** Fast, scalable batch-effect integration
* **SingleR / celldex:** Automated reference-based cell-type annotation
* **SingleCellExperiment:** Interoperability conversion
* **ggplot2 / ggpubr:** Data visualization & figure layout

---

## Results Summary
Before integration, cells clustered primarily by experimental treatment condition (`CTRL` vs `STIM`). Following Harmony integration, cell clusters aligned strictly by biological cell type (e.g., T cells, B cells, Monocytes) while preserving condition-specific expression shifts within each aligned cluster.
