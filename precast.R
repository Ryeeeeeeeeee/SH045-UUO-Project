```r
# Set the working directory and load required libraries
getwd()
setwd('./rds')
library(PRECAST)
library(Seurat)
library(dplyr)
library(ggplot2)

# Load the Seurat objects
UUOVH2 <- readRDS('./rds/002_UUO_Veh.rds')
UUOIH4 <- readRDS('./rds/003_UUO_SH045.rds')

# Combine the two datasets into a list
bc2 <- list(UUOVH2, UUOIH4)
head(bc2[[1]])

# Get the gene-by-spot read count matrices
countList <- lapply(bc2, function(x) x[["RNA"]]@counts)
M <- length(countList)

# Get the metadata of each spot for each data batch
metadataList <- lapply(bc2, function(x) x@meta.data)

# Check if "row" and "col" columns are present in metadata and display their head
for (r in 1:M) {
  meta_data <- metadataList[[r]]
  all(c("row", "col") %in% colnames(meta_data))  # Check if names are correct
  head(meta_data[, c("row", "col")])
}

# Ensure the row names of metadata in metaList are the same as the column names of count matrix in countList
for (r in 1:M) {
  row.names(metadataList[[r]]) <- colnames(countList[[r]])
}

# Create the Seurat list object
seuList <- list()
for (r in 1:M) {
  seuList[[r]] <- CreateSeuratObject(counts = countList[[r]], meta.data = metadataList[[r]], project = "AgingHeartPRECAST")
}

# Update the bc2 list with the Seurat objects and remove the temporary list
bc2 <- seuList
rm(seuList)

# Display the head of metadata for row and col
head(meta_data[, c("row", "col")])

# Create PRECASTObject
set.seed(2022)
PRECASTObj <- CreatePRECASTObject(bc2, project = "BC2", gene.number = 2000, selectGenesMethod = "SPARK-X",
                                  premin.spots = 20, premin.features = 20, postmin.spots = 1, postmin.features = 10)

# Check the Seurat list in PRECASTObj
PRECASTObj@seulist
PRECASTObj@seuList  # Check if seuList is null

# Add adjacency matrix list for a PRECASTObj object to prepare for PRECAST model fitting
PRECASTObj <- AddAdjList(PRECASTObj, platform = "Visium")

# Add a model setting in advance for a PRECASTObj object with verbose output
PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = FALSE, verbose = TRUE, int.model = NULL)

# Fit the PRECAST model given K
PRECASTObj <- PRECAST(PRECASTObj, K = 10)

# Backup the fitting results in resList
resList <- PRECASTObj@resList

# Select the best model
PRECASTObj <- SelectModel(PRECASTObj)

# Integrate spatial data
seuInt <- IntegrateSpaData(PRECASTObj, species = "Mouse")
seuInt

# Choose colors for clusters and plot spatial clusters
cols_cluster <- chooseColors(palettes_name = "Classic 20", n_colors = 13, plot_colors = TRUE)
p12 <- SpaPlot(seuInt, item = "cluster", batch = NULL, point_size = 4.5, cols = cols_cluster, combine = TRUE,
               nrow.legend = 7)
p12

# Plot spatial clusters without combining
pList <- SpaPlot(seuInt, item = "cluster", batch = NULL, point_size = 3, cols = cols_cluster, combine = FALSE,
                 nrow.legend = 7)
drawFigs(pList, layout.dim = c(2, 4), common.legend = TRUE, legend.position = "right", align = "hv")

# Add UMAP and t-SNE dimensional reductions
seuInt <- AddUMAP(seuInt)
seuInt <- AddTSNE(seuInt, n_comp = 2)

# Plot UMAP and t-SNE plots with clusters, batch, and group information
p1 <- dimPlot(seuInt, item = "cluster", point_size = 0.5, font_family = "serif", cols = cols_cluster,
              border_col = "gray10", nrow.legend = 14, legend_pos = "right")  
p2 <- dimPlot(seuInt, item = "batch", point_size = 0.5, font_family = "serif", legend_pos = "right")
p3 <- dimPlot(seuInt, item = "group", point_size = 0.5, font_family = "serif", legend_pos = "right")
drawFigs(list(p1, p2), layout.dim = c(1, 2), legend.position = "right", align = "hv")

# Find differentially expressed genes for all clusters
dat_deg <- FindAllMarkers(seuInt)

# Select top 10 markers for each cluster
n <- 10
dat_deg %>%
  group_by(cluster) %>%
  top_n(n = n, wt = avg_log2FC) -> top10

# Scale the data and downsample to 400 cells
seuInt <- ScaleData(seuInt)
seus <- subset(seuInt, downsample = 400)
color_id <- as.numeric(levels(Idents(seus)))

# Create a heatmap of the top markers
p1 <- doHeatmap(seus, features = top10$gene, cell_label = "Domain", grp_label = FALSE, grp_color = cols_cluster[color_id],
                pt_size = 6, slot = "scale.data") + 
  theme(legend.text = element_text(size = 10), legend.title = element_text(size = 13, face = "bold"), 
        axis.text.y = element_text(size = 5, face = "italic", family = "serif"))
p1



