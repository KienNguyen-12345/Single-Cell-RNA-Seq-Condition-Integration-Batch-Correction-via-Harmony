#Integrate single-cell RNA-Seq data in R using Harmony | Harmony in Seurat Workflow

setwd('/Users/nguyenkien/Desktop/Practice R/Lesson 5/Lesson 5')

#Harmony is a batch correction method:
# Work directly on low-dimensional linear embeddings (e.g. PCA coordinates)

#Study Design:
# Goal of study: to assess cell-type-specific changes in gene expression with interferon (IFN)-beta treatment using data from SeuratData
# Desgin: Peripheral blood mononuclear cells (PBMCs) from 8 lupus patients were split into simulated and control group
# & the simulated group was treated with IFN beta.

# Goal of our analysis: to integrate data by conditions, to overlay cells taht are similar in both conditions
devtools::install_github("satijalab/seurat-data")
n#Load libraries:
library(harmony)
library(Seurat)
library(SeuratData)
library(tidyverse)
library(ggplot2)
library(ggpubr)

#View available data from SeuratData package:
AvailableData()
#Install available data from SeuratData package:
InstallData("ifnb")

#Load downloaded dataset & Updatefor Seurat V5:

LoadData("ifnb")
ifnb <- UpdateSeuratObject(ifnb)
#QC and filtering:
#QC:
#1. Mitochondrial percent:
ifnb$percent.mt <- PercentageFeatureSet(ifnb, pattern = "^MT-")
view(ifnb@meta.data)

#2. Feature scatter:
FeatureScatter(ifnb, 
              feature1 = "nCount_RNA",
              feature2 = "nFeature_RNA",
              log = T) + geom_smooth(method = 'lm')
#Result: inconclusive
#3. Vln Plot:
VlnPlot(ifnb,
        features = c('nCount_RNA', 'nFeature_RNA'), log=T)
#3a. Box Plot:

ggboxplot(ifnb@meta.data, 
          x = "seurat_annotations", 
          y = "nCount_RNA", 
          fill = "seurat_annotations",
          x.text.angle = 45)
#4: Filtering:

ifnb_filtered <- subset(ifnb,
                        subset = nCount_RNA > 800)
view(ifnb_filtered@meta.data)
VlnPlot(ifnb_filtered,
        features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), log=T)

#Normal standard workflow:
#Normalize data:
ifnb_filtered <- NormalizeData(ifnb_filtered)
#Find variable feature:
ifnb_filtered <- FindVariableFeatures(ifnb_filtered)
#Scale data:
ifnb_filtered <- ScaleData(ifnb_filtered)

#RunPCA and find neighbors + clusters:
ifnb_filtered <- RunPCA(ifnb_filtered, features = VariableFeatures(ifnb_filtered))
ElbowPlot(ifnb_filtered)

ifnb_filtered <- FindNeighbors(ifnb_filtered, dims = 1:20)
ifnb_filtered <- FindClusters(ifnb_filtered)
view(ifnb_filtered@meta.data)

#Run UMAP:
ifnb_filtered <- RunUMAP(ifnb_filtered, dims = 1:20, reduction = 'pca')
#visualize:

before_1 <- DimPlot(ifnb_filtered, reduction = 'umap', group.by = ('seurat_annotations')) + ggtitle('Before harmony')
before <- DimPlot(ifnb_filtered, reduction = 'umap', group.by = ('stim'))+ ggtitle('Before harmony')
before

#Add cell types for labels:
#load libraries:
library(SingleCellExperiment)
library(celldex)
library(SingleR)
#1: download human cell reference database:
reference_database <- celldex::HumanPrimaryCellAtlasData()
#2: convert data from Seurat to SingleCellExperiment:
ifnb_SingleR <- as.SingleCellExperiment(ifnb_filtered)
#3: Predict cell type:
predict_cell_type <- SingleR(test = ifnb_SingleR,
                             ref = reference_database,
                             labels = reference_database$label.main)
#4. Add this variable into Seurat:
ifnb_filtered$cell_labels <- predict_cell_type$labels

#Visualize:
plot_with_cell_labels <- DimPlot(ifnb_filtered, 
                            reduction = 'umap',
                            group.by = c('stim','cell_labels'),
                            label = T,
                            repel = T)
plot_with_cell_labels


# Run Harmony:
ifnb.harmony <- ifnb_filtered %>%
  RunHarmony(group.by.vars = 'stim', plot_convergence = F)

#RunUMAP and clustering using harmony instead of PCA:
ifnb.harmony <- ifnb.harmony%>%
  RunUMAP(reduction = 'harmony', dims= 1:20)%>%
  FindNeighbors(reduction = 'harmony', dims = 1:20)%>%
  FindClusters(resolution = 0.5)

#Visualize:
after <- DimPlot(ifnb.harmony,
        reduction = 'umap',
        group.by = c('stim')) + ggtitle('After harmony')
after
after1 <- DimPlot(ifnb.harmony,
                  reduction = 'umap',
                  group.by = 'seurat_annotations') + ggtitle('After harmony')

before|after

before_1|after1
