library(dplyr)
library(Seurat)
library(patchwork)


setwd("L:/Lab-Wiestner/Jonathan/scDNA Sequencing/Tutorials/data/pbmc3k_filtered_gene_bc_matrices.tar/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19")

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "L:/Lab-Wiestner/Jonathan/scDNA Sequencing/Tutorials/data/pbmc3k_filtered_gene_bc_matrices.tar/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc