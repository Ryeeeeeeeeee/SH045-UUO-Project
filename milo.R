```r
# Install and load necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Clear the workspace and run garbage collection
rm(list = ls())
gc()

# Install Bioconductor packages
BiocManager::install("miloR")
BiocManager::install("scran")
BiocManager::install("MouseGastrulationData")

# Load required libraries
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(MouseGastrulationData)
library(Seurat)

# Load example data from miloR package
data("sim_trajectory", package = "miloR")

# Load the merged Seurat object
all_merge <- readRDS('./rds/UUOVH_vs_UUOIH_merge_DF_SCTharmony_celldenifition_allsubsets.rds')

# Convert Seurat object to SingleCellExperiment object
embryo_data <- as.SingleCellExperiment(all_merge)

# Initialize a Milo object with the SingleCellExperiment data
embryo_milo <- Milo(embryo_data)
embryo_milo

# Build the k-nearest neighbor graph
embryo_milo <- buildGraph(embryo_milo, k = 50, d = 30)

# Create neighborhoods of cells
embryo_milo <- makeNhoods(embryo_milo, prop = 0.1, k = 50, d = 30, refined = FALSE, reduced_dims = "pca.corrected")

# Plot histogram of neighborhood sizes
plotNhoodSizeHist(embryo_milo)

# Count cells in each neighborhood
embryo_milo <- countCells(embryo_milo, meta.data = as.data.frame(colData(embryo_milo)), sample = "orig.ident")
head(nhoodCounts(embryo_milo))

# Prepare experimental design dataframe
embryo_design <- data.frame(colData(embryo_milo))[, c("orig.ident", "group")]
embryo_design$orig.ident <- as.factor(embryo_design$orig.ident) 
embryo_design <- distinct(embryo_design)
rownames(embryo_design) <- embryo_design$orig.ident

# Calculate neighborhood distances
embryo_milo <- calcNhoodDistance(embryo_milo, d = 30)

# Perform differential abundance testing
da_results <- testNhoods(embryo_milo, design = ~ group, design.df = embryo_design)
head(da_results)

# Plot results of differential abundance testing
da_results %>%
  arrange(SpatialFDR) %>%
  head()
ggplot(da_results, aes(PValue)) + geom_histogram(bins = 50)
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) # Mark significance threshold (10% FDR)

# Build neighborhood graph
embryo_milo <- buildNhoodGraph(embryo_milo)

# Plot single-cell UMAP
umap_pl <- plotReducedDim(embryo_milo, dimred = "UMAP", colour_by = "group", text_by = "cell_type", 
                          text_size = 3, point_size = 0.5) +
  guides(fill = "none")

# Plot neighborhood graph
nh_graph_pl <- plotNhoodGraphDA(embryo_milo, da_results, layout = "UMAP", alpha = 0.1) 

# Combine UMAP and neighborhood graph plots
umap_pl + nh_graph_pl +
  plot_layout(guides = "collect")

## Extract SingleCellExperiment object from example data
traj_sce <- sim_trajectory[['SCE']]

## Extract sample metadata for testing
traj_meta <- sim_trajectory[["meta"]]

## Add metadata to colData slot of SingleCellExperiment
colData(traj_sce) <- DataFrame(traj_meta)
logcounts(traj_sce) <- log(counts(traj_sce) + 1)

# Perform PCA on the SingleCellExperiment object
traj_sce <- runPCA(traj_sce, ncomponents = 30)

# Perform UMAP on the SingleCellExperiment object
traj_sce <- runUMAP(traj_sce)

# Plot UMAP
plotUMAP(traj_sce)

# Initialize a Milo object with the trajectory data
traj_milo <- Milo(traj_sce)
reducedDim(traj_milo, "UMAP") <- reducedDim(traj_sce, "UMAP")

# Build k-nearest neighbor graph for trajectory data
traj_milo <- buildGraph(traj_milo, k = 10, d = 30)

# Create neighborhoods of cells for trajectory data
traj_milo <- makeNhoods(traj_milo, prop = 0.1, k = 10, d = 30, refined = TRUE)
plotNhoodSizeHist(traj_milo)

# Count cells in each neighborhood for trajectory data
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), samples = "Sample")
head(nhoodCounts(traj_milo))

# Prepare experimental design dataframe for trajectory data
traj_design <- data.frame(colData(traj_milo))[, c("Sample", "Condition")]
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$Sample

# Reorder rownames to match columns of neighborhood counts
traj_design <- traj_design[colnames(nhoodCounts(traj_milo)), , drop = FALSE]

# Calculate neighborhood distances for trajectory data
traj_milo <- calcNhoodDistance(traj_milo, d = 30)

# Perform differential abundance testing for trajectory data
da_results <- testNhoods(traj_milo, design = ~ Condition, design.df = traj_design)
da_results %>%
  arrange(-SpatialFDR) %>%
  head()

# Build neighborhood graph for trajectory data
traj_milo <- buildNhoodGraph(traj_milo)

# Plot UMAP and neighborhood graph for trajectory data
plotUMAP(traj_milo) + plotNhoodGraphDA(traj_milo, da_results, alpha = 0.05) +
  plot_layout(guides = "collect")
```