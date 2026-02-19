# PBMC3K Single-Cell RNA-seq Analysis
# Author: Adekunle Ajiboye
# Date: 02/18/2026

# Load libraries
library(Seurat)

# Load data
pbmc.data <- Read10X(data.dir = ".")

pbmc <- CreateSeuratObject(
  counts = pbmc.data,
  project = "pbmc3k",
  min.cells = 3,
  min.features = 200
)

pbmc

# Quality Control
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

VlnPlot(pbmc,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)

pbmc <- subset(pbmc,
               subset = nFeature_RNA > 200 &
                 nFeature_RNA < 2500 &
                 percent.mt < 5)

# Normalization
pbmc <- NormalizeData(pbmc)

# Find variable genes
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
length(VariableFeatures(pbmc))

# Scaling
pbmc <- ScaleData(pbmc)

# PCA
pbmc <- RunPCA(pbmc)
ElbowPlot(pbmc)

# Clustering
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, label = TRUE)

# Plot
markers <- FindAllMarkers(pbmc, only.pos = TRUE)
head(markers)

FeaturePlot(pbmc, features = c("CD3D", "MS4A1", "LYZ"))
