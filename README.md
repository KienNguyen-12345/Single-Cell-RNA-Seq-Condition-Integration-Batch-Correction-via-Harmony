

  <h1>Single-Cell RNA-Seq Integration and Differential Expression with Harmony &amp; Seurat</h1>
  <p>A reproducible R workflow demonstrating batch effect correction, cell type annotation, and condition-specific differential expression analysis using <strong>Seurat (v5)</strong> and <strong>Harmony</strong>[cite: 1].</p>

  <hr>

  <h2>Overview</h2>
  <ul>
    <li><strong>Dataset:</strong> Peripheral blood mononuclear cells (PBMCs) from systemic lupus erythematosus (SLE) patients split into interferon-beta stimulated (<code>STIM</code>) and unstimulated control (<code>CTRL</code>) conditions (<code>SeuratData::ifnb</code>).</li>
    <li><strong>Primary Objective:</strong> Integrate multi-condition scRNA-seq profiles across experimental batches to align shared cell states and identify condition-dependent transcriptional perturbations[cite: 1].</li>
    <li><strong>Key Analytical Approaches:</strong>
      <ul>
        <li>Low-dimensional embedding integration via <strong>Harmony</strong>.</li>
        <li>Automated reference-based cell type annotation using <strong>SingleR</strong> and <code>celldex</code>.</li>
        <li>Identification of cell-state conserved markers alongside treatment-responsive marker discovery using <code>DESeq2</code> and Seurat testing suites.</li>
      </ul>
    </li>
  </ul>

  <hr>

  <h2>Workflow Architecture</h2>
  <pre><code>Load IFNB Dataset
  │
  ▼
QC &amp; Filtering: nCount_RNA &gt; 800, MT%
  │
  ▼
Standard Seurat Workflow: Normalize, Scale, PCA
  ├────────────────────────┬────────────────────────┐
  ▼                        ▼                        ▼
Unintegrated DimPlot    Automated Annotation    RunHarmony: Correct
(Separation by 'stim')  via SingleR + HPCA      Embeddings by 'stim'
                                                    │
                                                    ▼
                                                Harmony Neighbors,
                                                Clusters &amp; UMAP
                                                    │
                                                    ▼
                                                Conserved Markers &amp;
                                                Differential Expression
                                                (STIM vs CTRL)</code></pre>

  <hr>

  <h2>Pipeline Steps</h2>

  <h3>1. Quality Control &amp; Preprocessing</h3>
  <ul>
    <li>Calculated mitochondrial gene percentage (<code>PercentageFeatureSet(pattern = "^MT-")</code>).</li>
    <li>Inspected feature-count distributions using <code>FeatureScatter</code>, <code>VlnPlot</code>, and metadata boxplots.</li>
    <li>Filtered low-quality cells retaining barcodes with <code>nCount_RNA &gt; 800</code>.</li>
  </ul>

  <h3>2. Standard Processing &amp; Unintegrated Baseline</h3>
  <ul>
    <li>Log-normalized gene expression (<code>NormalizeData</code>), selected top highly variable genes (<code>FindVariableFeatures</code>), and scaled variance (<code>ScaleData</code>).</li>
    <li>Computed PCA components and evaluated dimensionality with <code>ElbowPlot</code>.</li>
    <li>Generated baseline UMAP projections confirming significant batch/condition separation driven by interferon stimulation (<code>stim</code>).</li>
  </ul>

  <h3>3. Reference-Based Annotation</h3>
  <ul>
    <li>Converted the Seurat object to <code>SingleCellExperiment</code> format.</li>
    <li>Queried the <code>HumanPrimaryCellAtlasData</code> reference from <code>celldex</code> via <code>SingleR</code> to project automated cell identity labels onto single cells.</li>
  </ul>

  <h3>4. Batch Correction via Harmony</h3>
  <ul>
    <li>Applied <code>RunHarmony(group.by.vars = 'stim')</code> on PCA embeddings to align cells across stimulation conditions while preserving biological variance.</li>
    <li>Computed neighborhood graphs (<code>FindNeighbors</code>), clustered cells at <code>resolution = 0.5</code>, and projected onto Harmony-corrected UMAP space[cite: 1].</li>
  </ul>

  <h3>5. Marker Identification &amp; Cell Type Resolution</h3>
  <ul>
    <li><strong>All Cluster Markers:</strong> Identified positive enriched genes with <code>FindAllMarkers(test.use = 'DESeq2', logfc.threshold = 0.25, min.pct = 0.1)</code>[cite: 1].</li>
    <li><strong>Conserved Markers:</strong> Evaluated condition-invariant identity markers across <code>stim</code> groups using <code>FindConservedMarkers()</code> (identifying canonical CD4 T cell markers such as <em>HLA-DRA</em> and <em>HLA-DRB1</em>)[cite: 1].</li>
    <li>Renamed cluster identities to corresponding lineages (e.g., Cluster 3 annotated as <code>CD4 Naive T</code>)[cite: 1].</li>
  </ul>

  <h3>6. Condition-Specific Differential Expression</h3>
  <ul>
    <li>Generated composite metadata identities (<code>celltype.condition</code>, e.g., <code>CD4 Naive T_STIM</code> vs. <code>CD4 Naive T_CTRL</code>)[cite: 1].</li>
    <li>Executed <code>FindMarkers()</code> to isolate interferon-beta responsive genes (e.g., induction of <em>IFIT3</em>)[cite: 1].</li>
    <li>Contrasted conserved vs. induced markers across conditions using split-view <code>FeaturePlot</code> panels[cite: 1].</li>
  </ul>

  <hr>

  <h2>Dependencies</h2>
  <p>Ensure the following R packages are installed before running <code>Practice 5.R</code>[cite: 1]:</p>
  <pre><code># Bioconductor &amp; CRAN packages
install.packages(c("tidyverse", "ggplot2", "ggpubr", "devtools"))
BiocManager::install(c("SingleCellExperiment", "celldex", "SingleR", "DESeq2"))

# Harmony and Seurat packages
install.packages("Seurat")
install.packages("harmony")
devtools::install_github("satijalab/seurat-data")</code></pre>

  <hr>

  <h2>Usage</h2>
  <ol>
    <li>Open <code>Practice 5.R</code> in RStudio or your preferred IDE[cite: 1].</li>
    <li>Adjust your working directory if needed (<code>setwd(...)</code>)[cite: 1].</li>
    <li>Run through sections sequentially to reproduce the integration embeddings, cluster assignments, and differential expression plots[cite: 1].</li>
  </ol>

</body>
</html>
