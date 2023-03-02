library(dplyr)
library(Seurat)
library(patchwork)


setwd("L:/Lab-Wiestner/Jonathan/scDNA Sequencing/Tutorials/data/pbmc3k_filtered_gene_bc_matrices.tar/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19")


# Setup Seurat Object -----------------------------------------------------
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "L:/Lab-Wiestner/Jonathan/scDNA Sequencing/Tutorials/data/pbmc3k_filtered_gene_bc_matrices.tar/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc


#### Standard pre-processing workflow

# 1) QC & Selecting Cells for Downstream ----------------------------------
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

head(pbmc@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# 2) Normalizing Data -------------------------------------------------------
#normalizes feature expression measurement by total expression, multiplies by scaling factor, and log transform results
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# 3) Feature Selection -------------------------------------------------------

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# 4) Scaling Data ------------------------------------------------------------
# linear transformation prior to dimensional reduction that makes mean expression 0 and variance 1
# ensures that highly-expressed genes do not dominate
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# 5) Perform Linear Dimensional Reduction (PCA) ---------------------------
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# 6) Determine Dimensionality ---------------------------------------------
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

JackStrawPlot(pbmc, dims = 1:15)

ElbowPlot(pbmc)

# 7) Clustering Cells -----------------------------------------------------
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

head(Idents(pbmc), 5)

# 8) Non-Linear Dimensional Reduction (UMAP/tSNE) -------------------------
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")


# 9) Finding Differentially Expressed Features (Clustered Biomarke --------

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

VlnPlot(pbmc, features = c("CCR7", "CD79A"))


