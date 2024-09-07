library(tidyverse)
library(Seurat)
library(Matrix)
```r
# Load required libraries
library(cowplot)
library(biomaRt)

# Load the Seurat object
scRNA <- readRDS('F:/game/kidney/reviewer_atlas/endo.rds')

# Read the marker genes file and filter for cluster 4
marker <- read.csv('F:/kidney/output/UUOVH_vs_UUOIH/Endothelial/test/Endothelial.markersSeurat.csv', row.names = 1)
marker <- marker[marker$cluster == 4,]
marker <- head(marker, 50)

# Extract mouse gene symbols
musGenes <- rownames(marker)

# Set up Ensembl BioMart for human and mouse
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

# Map mouse genes to human orthologs using Ensembl BioMart
genes <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                values = musGenes, 
                mart = mouse, 
                attributesL = c("hgnc_symbol"), 
                martL = human, uniqueRows = TRUE)

# Extract human gene symbols
gene <- genes$HGNC.symbol
gene <- list(gene)

# Add module score for endothelial markers to Seurat object
seurat <- AddModuleScore(object = scRNA,
                         features = gene,
                         ctrl = 5,
                         name = 'endo')

# Inspect column names of meta data
colnames(seurat@meta.data)

# Load ggplot2 library for visualization
library(ggplot2)

# Fetch UMAP coordinates and module scores from Seurat object
mydata <- FetchData(seurat, vars = c('umap_1', 'umap_2', 'endo1'))

# Create a scatter plot colored by module score
a <- ggplot(mydata, aes(x = umap_1, y = umap_2, colour = endo1)) + 
  geom_point(size = 0.5) + 
  scale_color_gradientn(values = seq(0, 1, 0.2), colours = c('grey', 'red'))

# Customize plot appearance and add theme
a + theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = 'black'))

# Create a new metadata column to classify cells based on module score
seurat[['stage']] <- ifelse(seurat@meta.data[,'endo1'] > 0.3, 'High', 'Low')

# Set the new classification as the active identity class
Idents(seurat) <- 'stage'

# Plot UMAP colored by the new classification
DimPlot(seurat)
