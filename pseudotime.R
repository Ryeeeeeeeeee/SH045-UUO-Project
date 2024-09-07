
# Clear environment and set working directory
rm(list = ls())
setwd('E:/output/UUOVH_vs_UUOIH/Fibroblasts/')

# Load necessary libraries and check package versions
library(monocle3)
package.version("monocle3")
library(Seurat)
library(SeuratObject)
library(tidyselect)
library(dplyr)
library(BiocGenerics)

# Read the Seurat object containing Fibroblasts data
Fibroblasts <- readRDS("E:/output/UUOVH_vs_UUOIH/Fibroblasts/Fibroblasts_Seurat.rds")

# Extract RNA assay data and cell metadata
data <- GetAssayData(Fibroblasts, assay = 'RNA', slot = 'counts')
cell_metadata <- Fibroblasts@meta.data

# Create gene annotation dataframe with gene short names
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

# Create a new CellDataSet object (cds) for Monocle 3
cds <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)

# Preprocess the CellDataSet
cds <- preprocess_cds(cds, num_dim = 50)

# Align cells by 'plate' and model residuals
cds <- align_cds(cds, alignment_group = "plate", residual_model_formula_str = '~Size_Factor')

# Reduce dimensionality using PCA
cds <- reduce_dimension(cds, preprocess_method = "PCA")

# Extract column names from colData
colnames(colData(cds))

# Plot cells colored by Seurat clusters without labeling groups by cluster
plot_cells(cds, label_groups_by_cluster = FALSE, color_cells_by = "seurat_clusters")

# Plot UMAP plot for cds
p1 <- plot_cells(cds, reduction_method = "UMAP", color_cells_by = "seurat_clusters") + ggtitle('cds.umap')

# Align int.embed with cds.embed and plot UMAP for int.embed
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(Fibroblasts, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method = "UMAP", color_cells_by = "seurat_clusters") + ggtitle('int.umap')

# Combine p1 and p2 plots side by side
p1 | p2

# Plot cells with specific genes, without labeling cell groups and trajectory graph
plot_cells(cds, genes = ciliated_genes, label_cell_groups = FALSE, show_trajectory_graph = FALSE)

# Cluster cells using resolution parameter
cds <- cluster_cells(cds, resolution = 0.00000000000001)

# Plot cells colored by partition
plot_cells(cds, color_cells_by = "partition")

# Learn cell trajectory graph
cds <- learn_graph(cds)

# Plot cells colored by Seurat clusters, with labels for cell groups, leaves, and branch points
plot_cells(cds, color_cells_by = "seurat_clusters", label_cell_groups = FALSE, label_leaves = TRUE, label_branch_points = TRUE, graph_label_size = 1.5)

# Order cells based on pseudotime
cds <- order_cells(cds)

# Plot cells colored by pseudotime, without labels for leaves and branch points
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 1.5)

# Plot cells colored by pseudotime on the left and by cell type on the right, with labeled cell groups
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 1.5) | 
  plot_cells(cds, color_cells_by = "cell_type", label_cell_groups = TRUE, label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 1.5)
