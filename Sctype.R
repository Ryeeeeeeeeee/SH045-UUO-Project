# Clear environment, perform garbage collection, and check working directory
rm(list = ls())
gc()
getwd()

# Load necessary libraries
library(tidyverse)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(readr)
library(cowplot)
library(openxlsx)
library(monocle3)
library(HGNChelper)
library(data.table)
library(magrittr)
library(ggpubr)
library(paletteer)
library(gplots)
library(ggsci)
library(stringr)

# Read Seurat object from file
seuratfile <- readRDS("F:/kidney/output/UUOVH_vs_UUOIH_vs_CTRL_merge_DF_SCTharmony_sctype.rds")

# Display names of meta data columns in Seurat object
names(seuratfile@meta.data)

# Load function for preparing gene sets from GitHub
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")

# Load function for scoring cell types from GitHub
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# Define URL for database file and tissue of interest
db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
tissue <- c("Kidney")

# Prepare gene sets for the specified tissue
gs_list <- gene_sets_prepare(db_, tissue)

# Score cell types using prepared gene sets
es.max <- sctype_score(scRNAseqData = seuratfile[["SCT"]]@scale.data, scaled = TRUE, 
                       gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# Merge scores by cluster and summarize top cell types
cL_results <- do.call("rbind", lapply(unique(seuratfile@meta.data$seurat_clusters), function(cl){
  es.max.cl <- sort(rowSums(es.max[, rownames(seuratfile@meta.data[seuratfile@meta.data$seurat_clusters == cl, ])]), decreasing = TRUE)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seuratfile@meta.data$seurat_clusters == cl)), 10)
}))

# Select top cell type per cluster based on scores
sctype_scores <- cL_results %>% 
  group_by(cluster) %>% 
  top_n(n = 1, wt = scores)

# Assign 'Unknown' to cell types with scores less than 25% of the total cells
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"

# Print and write the summarized results to a CSV file
print(sctype_scores[, 1:3])
write.csv(sctype_scores, "scores_all_inh.csv", row.names = FALSE, quote = FALSE)

# Assign custom classification to Seurat object based on top cell types per cluster
seuratfile@meta.data$customclassif <- ""
for (j in unique(sctype_scores$cluster)) {
  cl_type <- sctype_scores[sctype_scores$cluster == j, ]
  seuratfile@meta.data$customclassif[seuratfile@meta.data$seurat_clusters == j] <- as.character(cl_type$type[1])
}

# Create a 'group' factor variable from 'orig.ident'
seuratfile$group <- as.factor(gsub("[[:digit:]]", "", seuratfile$orig.ident))

# Plot UMAP with labels and grouping by 'group'
DimPlot(seuratfile, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'group')

# Rename and recode custom classifications for clarity
seuratfile$customclassif <- gsub('Hematopoietic cells|Juxtaglomerular cells|Immune cells', 'Inflammatory cells', seuratfile$customclassif)
seuratfile$customclassif <- gsub('α-intercalated', 'Intercalated', seuratfile$customclassif)
seuratfile$customclassif <- gsub('Stromal cells', 'Fibroblasts', seuratfile$customclassif)
unique(seuratfile$customclassif)

# Recode custom classifications directly in 'sctype_scores'
sctype_scores$type <- gsub('Hematopoietic cells|Juxtaglomerular cells|Immune cells', 'Inflammatory cells', sctype_scores$type)
sctype_scores$type <- gsub('α-intercalated', 'Intercalated', sctype_scores$type)
sctype_scores$type <- gsub('Stromal cells', 'Fibroblasts', sctype_scores$type)

# Arrange and identify new cluster IDs
sctype_scores <- sctype_scores %>% arrange(cluster)
Idents(seuratfile) <- ('customclassif')

# Rename cluster identities using new cluster IDs
new.cluster.ids <- as.vector(sctype_scores$type)
names(new.cluster.ids) <- levels(seuratfile)
seuratfile <- RenameIdents(seuratfile, new.cluster.ids)

# Plot UMAP with labels colored by 'customclassif'
DimPlot(seuratfile, label = TRUE)
