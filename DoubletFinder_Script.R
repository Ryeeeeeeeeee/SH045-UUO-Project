```r
# Clear the workspace
rm(list = ls())

# Load required libraries
library(dplyr)
library(Seurat)
library(patchwork)
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
library(ggplot2)

# Set the working directory
setwd("E:/output/UUOIH1")

# Read the Seurat object
UUOIH1 <- readRDS("E:/game/output/UUOIH1/UUOIH1_Seurat.rds")

###### DoubletFinder Analysis #######

# Step 1: Parametric Sweep for DoubletFinder
# Perform a parametric sweep using the first 15 principal components (PCs)
sweep.res.list_UUOIH1 <- paramSweep_v3(UUOIH1, PCs = 1:15, sct = FALSE)

# Summarize the results of the parametric sweep
sweep.stats_UUOIH1 <- summarizeSweep(sweep.res.list_UUOIH1, GT = FALSE)

# Find the optimal pK value
bcmvn <- find.pK(sweep.stats_UUOIH1)

# Step 2: Visualize the pK value
# Create a PNG file to save the plot
png(file = "E:/output/UUOIH1/UUOIH1_DoubletFinder_pK.png", width = 1040, height = 663)

# Extract pK and BCmetric values
pK = as.numeric(as.character(bcmvn$pK))
BCmetric = bcmvn$BCmetric

# Choose the pK value with the maximum BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

# Plot the BCmetric values and mark the chosen pK value
par(mar = c(5, 4, 4, 8) + 1, cex.main = 1.2, font.main = 2)
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "blue", lty = 1)
abline(v = pK_choose, lwd = 2, col = 'red', lty = 2)
title("The BCmvn Distributions")
text(pK_choose, max(BCmetric), as.character(pK_choose), pos = 4, col = "red")

# Close the PNG device
dev.off()

# Print the chosen pK value
pK_choose

###### Homotypic Doublet Proportion Estimate ########
# Calculate the obtained_vs_expected ratio for the given dataset
# The values are based on the number of cells obtained and expected for different samples
obtained_vs_expected <- 0.0390

# Estimate the homotypic doublet proportion using the Seurat clusters
homotypic.prop <- modelHomotypic(UUOIH1@meta.data$seurat_clusters)

# Calculate the expected number of doublets
nExp_poi <- round(obtained_vs_expected * nrow(UUOIH1@meta.data))

# Adjust the expected number of doublets by the homotypic proportion
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

###### DoubletFinder Analysis #####
# Run DoubletFinder with the estimated parameters
UUOIH1 <- doubletFinder_v3(seu = UUOIH1, PCs = 1:15, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE)

# Save the updated Seurat object
saveRDS(UUOIH1, "E:/output/UUOIH1/UUOIH1_singlet_DoubletFinder.rds")
