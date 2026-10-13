Single-Cell RNA-Seq Integration and Differential Expression with Harmony & Seurat

A reproducible R workflow demonstrating batch effect correction, cell type annotation, and condition-specific differential expression analysis using Seurat (v5) and Harmony.

Overview

Dataset: Peripheral blood mononuclear cells (PBMCs) from systemic lupus erythematosus (SLE) patients split into interferon-beta stimulated (STIM) and unstimulated control (CTRL) conditions (SeuratData::ifnb).

Primary Objective: Integrate multi-condition scRNA-seq profiles across experimental batches to align shared cell states and identify condition-dependent transcriptional perturbations.

Key Analytical Approaches:

Low-dimensional embedding integration via Harmony.

Automated reference-based cell type annotation using SingleR and celldex.

Identification of cell-state conserved markers alongside treatment-responsive marker discovery using DESeq2 and Seurat testing suites.

Workflow Architecture

flowchart TD
    A[Load IFNB Dataset] --> B[QC & Filtering: nCount_RNA > 800, MT%]
    B --> C[Standard Seurat Workflow: Normalize, Scale, PCA]
    C --> D[Unintegrated DimPlot: Separation by 'stim']
    C --> E[Automated Annotation via SingleR + HPCA]
    C --> F[RunHarmony: Correct Embeddings by 'stim']
    F --> G[Harmony-based Neighbors, Clusters & UMAP]
    G --> H[Conserved Marker Discovery: FindConservedMarkers]
    G --> I[Condition Split: celltype.condition]
    I --> J[Differential Expression: FindMarkers STIM vs CTRL]


Pipeline Steps

1. Quality Control & Preprocessing

Calculated mitochondrial gene percentage (PercentageFeatureSet(pattern = "^MT-")).

Inspected feature-count distributions using FeatureScatter, VlnPlot, and metadata boxplots.

Filtered low-quality cells retaining barcodes with nCount_RNA > 800.

2. Standard Processing & Unintegrated Baseline

Log-normalized gene expression (NormalizeData), selected top highly variable genes (FindVariableFeatures), and scaled variance (ScaleData).

Computed PCA components and evaluated dimensionality with ElbowPlot.

Generated baseline UMAP projections confirming significant batch/condition separation driven by interferon stimulation (stim).

3. Reference-Based Annotation

Converted the Seurat object to SingleCellExperiment format.

Queried the HumanPrimaryCellAtlasData reference from celldex via SingleR to project automated cell identity labels onto single cells.

4. Batch Correction via Harmony

Applied RunHarmony(group.by.vars = 'stim') on PCA embeddings to align cells across stimulation conditions while preserving biological variance.

Computed neighborhood graphs (FindNeighbors), clustered cells at resolution = 0.5, and projected onto Harmony-corrected UMAP space.

5. Marker Identification & Cell Type Resolution

All Cluster Markers: Identified positive enriched genes with FindAllMarkers(test.use = 'DESeq2', logfc.threshold = 0.25, min.pct = 0.1).

Conserved Markers: Evaluated condition-invariant identity markers across stim groups using FindConservedMarkers() (identifying canonical CD4 T cell markers such as HLA-DRA and HLA-DRB1).

Renamed cluster identities to corresponding lineages (e.g., Cluster 3 annotated as CD4 Naive T).

6. Condition-Specific Differential Expression

Generated composite metadata identities (celltype.condition, e.g., CD4 Naive T_STIM vs. CD4 Naive T_CTRL).

Executed FindMarkers() to isolate interferon-beta responsive genes (e.g., induction of IFIT3).

Contrasted conserved vs. induced markers across conditions using split-view FeaturePlot panels.

Dependencies

Ensure the following R packages are installed before running Practice 5.R:

# Bioconductor & CRAN packages
install.packages(c("tidyverse", "ggplot2", "ggpubr", "devtools"))
BiocManager::install(c("SingleCellExperiment", "celldex", "SingleR", "DESeq2"))

# Harmony and Seurat packages
install.packages("Seurat")
install.packages("harmony")
devtools::install_github("satijalab/seurat-data")


Usage

Open Practice 5.R in RStudio or your preferred IDE.

Adjust your working directory if needed (setwd(...)).

Run through sections sequentially to reproduce the integration embeddings, cluster assignments, and differential expression plots.
