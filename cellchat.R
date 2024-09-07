```r
getwd()  # Get current working directory
setwd("E:/output/UUOIH/")  # Set working directory

rm(list = ls())  # Clear the environment
gc()  # Garbage collection to free up memory

# Suppress messages and install packages if they are not already installed
suppressMessages(if(!require(CellChat)) devtools::install_github("sqjin/CellChat"))
suppressMessages(if(!require(ggplot2)) install.packages('ggplot2'))
suppressMessages(if(!require(patchwork)) install.packages('patchwork'))
suppressMessages(if(!require(ggalluvial)) install.packages('ggalluvial'))
suppressMessages(if(!require(igraph)) install.packages('igraph'))
suppressMessages(if(!require(dplyr)) install.packages('dplyr'))

suppressMessages(options(stringsAsFactors = FALSE))
suppressMessages(options(future.globals.maxSize = 2*1024^3))  # Set max size for future global variables
suppressWarnings(suppressMessages(future::plan("multiprocess", workers = 12)))  # Set up parallel processing

library("CellChat")  # Load the CellChat library

# Load the Seurat object
seurat_object <- readRDS("E:/output/UUOIH/UUOIH_forcellchat_removeF2.rds")

# Extract normalized data matrix
data.input <- GetAssayData(seurat_object, assay = "RNA", slot = "data")

# Display column names of meta data
names(seurat_object@meta.data)

# Set cell identities
Idents(seurat_object) <- 'cell_type2'

# Create a dataframe of the cell labels
meta <- seurat_object@meta.data

# Create CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type2")

# Show factor levels of the cell labels
levels(cellchat@idents)

# Get the number of cells in each cell group
groupSize <- as.numeric(table(cellchat@idents))
groupSize

# Use the default CellChat mouse database
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

# Display the structure of the CellChat database interactions
dplyr::glimpse(CellChatDB$interaction)

# Assign the CellChat database to the CellChat object
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

# Subset data
cellchat <- subsetData(cellchat, features = NULL)

# Identify overexpressed genes
cellchat <- identifyOverExpressedGenes(cellchat)

# Identify overexpressed interactions
cellchat <- identifyOverExpressedInteractions(cellchat)

# Project data onto protein-protein interaction (PPI) network
cellchat <- projectData(cellchat, PPI.mouse)

# Compute communication probability
cellchat <- computeCommunProb(cellchat, raw.use = T)

# Filter communication
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Aggregate networks
cellchat <- aggregateNet(cellchat)

# Subset communication
df.net <- subsetCommunication(cellchat)
write.csv(df.net, '01.df.net.csv')

# Compute communication probability for pathways
cellchat <- computeCommunProbPathway(cellchat)

# Subset communication for pathways
df.netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netP, '01.df.netP.csv')

# Subset communication from specific sources
df.net.subset <- subsetCommunication(cellchat, sources.use = c('Endothelial_cell'))
write.csv(df.net.subset, 'Endothelial.df.net.csv')

# Display group sizes
groupSize

# Load patchwork library
library(patchwork)

# Plot interaction networks
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Interaction weights/strength")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Interaction weights/strength", sources.use = 'Endothelial cells')

# Plot interaction networks for each pathway
mat <- cellchat@net$weight
par(mfrow = c(3, 3), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# Display all pathway names
pathways.show <- df.netP$pathway_name

# Hierarchy plot
par(mfrow = c(1, 2), xpd = TRUE)
vertex.receiver <- seq(1, 5)
netVisual_aggregate(cellchat, signaling = pathways.show[3], vertex.receiver = vertex.receiver, layout = 'hierarchy')

# Circle plot
netVisual_aggregate(cellchat, signaling = pathways.show[3], layout = "circle")

# Create directory for circle plots
dir.create('E:/output/UUOIH/cellchat/circle')

# Save circle plots for all pathways
for (i in 1:length(pathways.show)) {
  pdf(paste0('E:/output/UUOIH/cellchat/circle/', pathways.show[i], '.pdf'), width = 7, height = 7)
  netVisual_aggregate(cellchat, signaling = pathways.show[i], layout = "circle")
  dev.off()
}

# Display cell identity levels
levels(cellchat@idents)

# Set vertex.receiver and create directory for pathway contribution plots
vertex.receiver <- seq(1, 4)
dir.create('pathwat.show/')

# Save pathway contribution plots
for (i in 1:length(pathways.show.all)) {
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename = paste0('pathwat.show/', pathways.show.all[i], "_L-R_contribution.pdf"), plot = gg, width = 10, height = 10, units = 'in', dpi = 300)
}

# Compute centrality of signaling roles
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Plot signaling role network
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 15, height = 6, font.size = 10)

# Select optimal number of patterns
selectK(cellchat, pattern = "outgoing")

# Save the CellChat object
saveRDS(cellchat, "E:/output/UUOIH/cellchat/UUOIH_cellchat.rds")
