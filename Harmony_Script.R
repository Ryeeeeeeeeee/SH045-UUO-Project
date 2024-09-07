```r
###### Harmony Integration and UMAP Visualization ######

# Clear the workspace
rm(list = ls())

# Get the current working directory
getwd()

# Set the working directory
setwd("E:/")

# Load required libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)
library(cowplot)

### Processing UUOIH Samples ###

# Load Seurat objects for UUOIH samples after DoubletFinder
UUOIH1 <- readRDS("E:/output/UUOIH1/UUOIH1_singlet_DoubletFinder.rds")
UUOIH2 <- readRDS("E:/output/UUOIH2/UUOIH2_singlet_DoubletFinder.rds")
UUOIH3 <- readRDS("E:/output/UUOIH3/UUOIH3_singlet_DoubletFinder.rds")

# Create lists of Seurat objects and their corresponding names
list_VRP1 <- c("UUOIH1", "UUOIH2", "UUOIH3")
list_VRP2 <- c(UUOIH1, UUOIH2, UUOIH3)

# Add Doublets column to metadata for each Seurat object
for (i in 1:length(list_VRP2)) {
  list_VRP2[[i]]$Doublets <- list_VRP2[[i]]@meta.data[, grep("DF.", colnames(list_VRP2[[i]]@meta.data))]
  assign(list_VRP1[i], list_VRP2[[i]])
}

# Merge UUOIH samples into a single Seurat object
UUOIH_merge <- merge(UUOIH1, y = c(UUOIH2, UUOIH3), project = "UUOIH_pipeline")
saveRDS(UUOIH_merge, "E:/output/UUOIH_merge.rds")

# Perform SCT normalization on the merged object
UUOIH_merge <- SCTransform(UUOIH_merge, verbose = FALSE)

# Run PCA on the SCT assay
UUOIH_merge <- RunPCA(object = UUOIH_merge, assay = "SCT")

# Run Harmony integration
UUOIH_merge <- UUOIH_merge %>%
  RunHarmony("orig.ident", assay.use = "SCT", plot_convergence = TRUE)

# Find neighbors and clusters, then run UMAP
UUOIH_merge <- FindNeighbors(UUOIH_merge, assay = "SCT", reduction = "harmony", dims = 1:20)
UUOIH_merge <- FindClusters(UUOIH_merge, resolution = 0.5)
UUOIH_merge <- RunUMAP(UUOIH_merge, assay = "SCT", reduction = "harmony", dims = 1:20)

# Save UMAP plots to PDF
pdf(file = "E:/output/UUOIH_UMAP.pdf")
DimPlot(UUOIH_merge, label = TRUE)
DimPlot(UUOIH_merge, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size = 0.5)
dev.off()

# Save the processed Seurat object
saveRDS(UUOIH_merge, "E:/output/UUOIH_merge_DF_SCTharmony.rds")

### Processing UUOVH Samples ###

# Load Seurat objects for UUOVH samples after DoubletFinder
UUOVH1 <- readRDS("E:/output/UUOVH1/UUOVH1_singlet_DoubletFinder.rds")
UUOVH2 <- readRDS("E:/output/UUOVH2/UUOVH2_singlet_DoubletFinder.rds")
UUOVH3 <- readRDS("E:/output/UUOVH3/UUOVH3_singlet_DoubletFinder.rds")

# Create lists of Seurat objects and their corresponding names
list_VRP1 <- c("UUOVH1", "UUOVH2", "UUOVH3")
list_VRP2 <- c(UUOVH1, UUOVH2, UUOVH3)

# Add Doublets column to metadata for each Seurat object
for (i in 1:length(list_VRP2)) {
  list_VRP2[[i]]$Doublets <- list_VRP2[[i]]@meta.data[, grep("DF.", colnames(list_VRP2[[i]]@meta.data))]
  assign(list_VRP1[i], list_VRP2[[i]])
}

# Merge UUOVH samples into a single Seurat object
UUOVH_merge <- merge(UUOVH1, y = c(UUOVH2, UUOVH3), project = "UUOVH_pipeline")
saveRDS(UUOVH_merge, "E:/output/UUOVH_merge.rds")

# Perform SCT normalization on the merged object
UUOVH_merge <- SCTransform(UUOVH_merge, verbose = FALSE)

# Run PCA on the SCT assay
UUOVH_merge <- RunPCA(object = UUOVH_merge, assay = "SCT")

# Run Harmony integration
UUOVH_merge <- UUOVH_merge %>%
  RunHarmony("orig.ident", assay.use = "SCT", plot_convergence = TRUE)

# Find neighbors and clusters, then run UMAP
UUOVH_merge <- FindNeighbors(UUOVH_merge, assay = "SCT", reduction = "harmony", dims = 1:20)
UUOVH_merge <- FindClusters(UUOVH_merge, resolution = 0.5)
UUOVH_merge <- RunUMAP(UUOVH_merge, assay = "SCT", reduction = "harmony", dims = 1:20)

# Save UMAP plots to PDF
pdf(file = "E:/output/UUOVH_UMAP.pdf")
DimPlot(UUOVH_merge, label = TRUE)
DimPlot(UUOVH_merge, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size = 0.5)
dev.off()

# Save the processed Seurat object
saveRDS(UUOVH_merge, "E:/output/UUOVH_merge_DF_SCTharmony.rds")

### Merging UUOVH and UUOIH Samples ###

# Create lists of all Seurat objects and their corresponding names
list_all1 <- c("UUOVH1", "UUOVH2", "UUOVH3", "UUOIH1", "UUOIH2", "UUOIH3")
list_all2 <- c(UUOVH1, UUOVH2, UUOVH3, UUOIH1, UUOIH2, UUOIH3)

# Add Doublets column to metadata for each Seurat object
for (i in 1:length(list_all2)) {
  list_all2[[i]]$Doublets <- list_all2[[i]]@meta.data[, grep("DF.", colnames(list_all2[[i]]@meta.data))]
  assign(list_all1[i], list_all2[[i]])
}

# Merge all samples into a single Seurat object
all_merge <- merge(UUOVH1, y = c(UUOVH2, UUOVH3, UUOIH1, UUOIH2, UUOIH3), project = "UUOVH_vs_UUOIH")

# Perform SCT normalization on the merged object
all_merge <- SCTransform(all_merge, verbose = FALSE)

# Run PCA on the SCT assay
all_merge <- RunPCA(object = all_merge, assay = "SCT")

# Run Harmony integration
all_merge <- all_merge %>%
  RunHarmony("orig.ident", assay.use = "SCT", plot_convergence = TRUE)

# Find neighbors and clusters, then run UMAP
all_merge <- FindNeighbors(all_merge, assay = "SCT", reduction = "harmony", dims = 1:30)
all_merge <- FindClusters(all_merge, resolution = 0.5)
all_merge <- RunUMAP(all_merge, assay = "SCT", reduction = "harmony", dims = 1:30)

# Display cluster identities
head(Idents(all_merge))

# Save UMAP plots to PDF
pdf(file = "E:/output/UUOVH_vs_UUOIH_UMAP.pdf")
DimPlot(all_merge, label = TRUE)
DimPlot(all_merge, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size = 0.5)
dev.off()

# Save the final processed Seurat object
saveRDS(all_merge, "E:/output/UUOVH_vs_UUOIH_merge_DF_SCTharmony.rds")




