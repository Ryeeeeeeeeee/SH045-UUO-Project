```r
rm(list = ls())  # Clear the environment
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
getwd()  # Get current working directory
setwd("E:/output/UUOIH1")  # Set working directory

# Load data with Read10X()
data <- Read10X(data.dir = "E:/game/UUOIH1_matrix_10X")

# Initialize the Seurat object with the raw (non-normalized data)
UUOIH1 <- CreateSeuratObject(counts = data, project = "UUOIH1", min.cells = 3, min.features = 200)

# Include percentage of mitochondrial genes expressed in your Seurat object metadata.
# This function recognizes gene names starting by "MT-" (if you are using mouse data, make sure to change it for "^mt-")
UUOIH1[["percent.mt"]] <- PercentageFeatureSet(UUOIH1, pattern = "^mt-")

# Visualize QC metrics as a violin plot
pdf(file = "E:/output/UUOIH1/UUOIH1_VlnPlot_FeatureScatterSeurat.pdf")  # File name

VlnPlot(UUOIH1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter to check number of counts vs number of genes or vs percent.mt
FeatureScatter(UUOIH1, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(UUOIH1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Closing the graphical device
dev.off()

# Subset by number of features and mitochondrial genes percentage
UUOIH1 <- subset(UUOIH1, subset = percent.mt < 10 & nFeature_RNA < 5000)

# Normalizing the data
UUOIH1 <- NormalizeData(UUOIH1, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features (feature selection)
UUOIH1 <- FindVariableFeatures(UUOIH1, selection.method = "vst", nfeatures = 1500)

# Scaling the data
all.genes <- rownames(UUOIH1)
UUOIH1 <- ScaleData(UUOIH1, features = all.genes)

write.csv(all.genes, "E:/output/UUOIH1/UUOIH1.all.genesSeurat.csv")

# Perform linear dimensional reduction
UUOIH1 <- RunPCA(UUOIH1, features = VariableFeatures(object = UUOIH1))
print(UUOIH1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(UUOIH1, dims = 1:2, reduction = "pca")
DimPlot(UUOIH1, reduction = "pca")
DimHeatmap(UUOIH1, dims = 1:15, cells = 500, balanced = TRUE)
UUOIH1 <- JackStraw(UUOIH1, num.replicate = 100)
UUOIH1 <- ScoreJackStraw(UUOIH1, dims = 1:20)
JackStrawPlot(UUOIH1, dims = 1:20)

# ElbowPlot
pdf(file = "E:/output/UUOIH1/UUOIH1_ElbowplotSeurat.pdf")  # File name

ElbowPlot(UUOIH1, ndims = 50)

# Closing the graphical device
dev.off()

# Calculate percentage of variability recapitulated by each component
# Get the total variance
eigValues = (UUOIH1[["pca"]]@stdev)^2  # EigenValues
varExplained2 = eigValues / as.numeric(slot(UUOIH1[["pca"]], "misc"))

pdf(file = "E:/output/UUOIH1/UUOIH1_VarExplained2Seurat.pdf")  # File name

plot(varExplained2)

# Closing the graphical device
dev.off()

# Cluster the cells
UUOIH1 <- FindNeighbors(UUOIH1, dims = 1:15)
UUOIH1 <- FindClusters(UUOIH1, resolution = 0.5)

# Run non-linear dimensional reduction (UMAP)
UUOIH1 <- RunUMAP(UUOIH1, dims = 1:15)

# Plot the UMAP
DimPlot(UUOIH1, reduction = "umap")
pdf(file = "E:/output/UUOIH1/UUOIH1_UMAPSeurat.pdf")  # File name

DimPlot(UUOIH1, reduction = "umap", label = TRUE)

# Closing the graphical device
dev.off()

# Quality check for the nFeatures and nGenes
pdf(file = "E:/output/UUOIH1/UUOIH1_violine clusterSeurat.pdf")  # File name

VlnPlot(UUOIH1, features = "nCount_RNA", idents = UUOIH1$seurat_clusters, pt.size = 0)
VlnPlot(UUOIH1, features = "nFeature_RNA", idents = UUOIH1$seurat_clusters, pt.size = 0)
VlnPlot(UUOIH1, features = "nCount_RNA", idents = UUOIH1$seurat_clusters)
VlnPlot(UUOIH1, features = "nFeature_RNA", idents = UUOIH1$seurat_clusters)

# Closing the graphical device
dev.off()

# Calculate the markers and store them in a new slot of @misc, which is empty by default and can be used to store different parameters
UUOIH1.markers <- FindAllMarkers(UUOIH1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
UUOIH1@misc$UUOIH11st.markers <- UUOIH1.markers

# Heatmap with top 10 markers for each cluster
pdf(file = "E:/output/UUOIH1/UUOIH1_Heatmaptop10genesSeurat.pdf")  # File name

top10 <- UUOIH1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(UUOIH1, features = top10$gene) + NoLegend()
DoHeatmap(UUOIH1, features = top10$gene, size = 1.5, group.bar.height = 0.01) + 
  theme(text = element_text(size = 8)) + theme(axis.text.y = element_text(size = 4))

# Closing the graphical device
dev.off()

# Get mean nFeatures and nCounts
p <- UUOIH1@meta.data %>% group_by(seurat_clusters) %>% 
  summarise(mean_nFeature_RNA = mean(nFeature_RNA), mean_nCount_RNA = mean(nCount_RNA)) %>% t()
p <- as.data.frame(p, row.names = c("", "mean_nFeatures", "mean_nCount"))

write.csv(p, "E:/output/UUOIH1/UUOIH1.mean_nFeat_nCountSeurat.csv")
write.csv(UUOIH1.markers, "E:/output/UUOIH1/UUOIH1.markersSeurat.csv")

# Save object to avoid following the whole process every time you need to access the data
saveRDS(UUOIH1, file = "E:/output/UUOIH1/UUOIH1_Seurat.rds")
```

