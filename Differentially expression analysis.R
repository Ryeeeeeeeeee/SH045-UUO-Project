```r
# Clear the workspace and set the working directory
rm(list = ls())
getwd()
setwd("E:/output/UUOVH_vs_UUOIH")

########################
### DESeq2 analysis for single cells
##########################

# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleCellExperiment")
library(DESeq2)
library(Seurat)
library(scran)
library(SingleCellExperiment)

#### Step 1
### Read the rds file
all_refined <- readRDS("E:/output/UUOVH_vs_UUOIH/UUOVH_vs_UUOIH_merge_DF_SCTharmony_celldenifition.rds")

### Convert the rds file to a matrix
all_refined1 = as.matrix(all_refined@assays$RNA@counts)

### Create a SingleCellExperiment object
sce <- SingleCellExperiment(all_refined1)

### Extract count data using assay() function
counts <- assay(sce)
counts(sce) <- counts

### Normalize the data using computeSumFactors from scran package
### Reference: https://bioconductor.statistik.tu-dortmund.de/packages/3.7/bioc/vignettes/scran/inst/doc/scran.html
sce <- computeSumFactors(sce)

### Log transform the data
logcounts(sce) <- log2(counts + 1)

### Optional: Save the SingleCellExperiment object
### saveRDS(sce, file = "./../data/all_resolution_0.7_sparsematrix_computesumfactors.rds")

gc()

### Extract the count matrix and log-transformed data
count_data <- assay(sce)
logcount_data <- assay(sce, "logcounts")

## -------------------- Creating meta data --------------------------

### Retrieve metadata from the original object
metadata <- all_refined@meta.data

gc()

### Modify group names
group <- metadata$group
group <- as.character(group)
group <- factor(group, levels = c("UUOIH", "UUOVH"))

### Create a data frame for coldata
coldata <- data.frame(cells = rownames(metadata), groups = group, cell_type = metadata$celltype)

gc()

## -------------------- Differentially expression analysis --------------------------
### Reference: https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html

### Loop through selected cell types for DESeq2 analysis
for(cc in c('Endothelial_cell','Fibroblasts')) { ### Specify clusters for analysis
  coldata <- data.frame(cells = rownames(metadata), groups = group, cell_type = metadata$seurat_clusters)
  gc()
  message(paste("Analysis starts for :", cc))
  coldata <- coldata[coldata$cell_type == cc,]
  
  ### Select count data for the specific cell type
  countdata <- assay(sce, "counts")
  countdata <- countdata[, colnames(countdata) %in% coldata$cells]
  
  ### Pre-filtering: Remove genes with low counts
  countdata <- countdata[rowSums(countdata) > 0,]
  
  message("Performing DESeqDataSetFromMatrix")
  ddsMat <- DESeqDataSetFromMatrix(countData = round(countdata + 1), 
                                   colData = coldata, 
                                   design = ~ groups)
  
  gc()
  
  ### Estimate size factors and dispersions
  dds <- estimateSizeFactors(ddsMat)
  dds <- estimateDispersionsGeneEst(dds)
  dispersions(dds) <- mcols(dds)$dispGeneEst
  
  message("Performing DESeq")
  dds <- DESeq(dds, fitType = "parametric")
  
  ### Create a directory for output
  celltype <- paste("cell_type", cc, sep = "_")
  dir.create('Differentially expression analysis UUOVH vs UUOIH/', celltype)
  outputdir <- paste('Differentially expression analysis UUOVH vs UUOIH/', celltype)
  gc()
  
  ### Output results of Wald test for contrast
  contrast <- c("groups", "UUOIH", "UUOVH")
  
  message("Performing results")
  res <- results(dds, contrast = contrast, alpha = 0.05)
  
  ### Write the results to a CSV file
  write.csv(as.data.frame(res), file = paste('Differentially expression analysis UUOVH vs UUOIH/', "/deseq2_cluster_rawcount_", cc, ".csv", sep = ""),
            row.names = TRUE)
  
  ### Plot MA and save as PDF
  pdf(paste('Differentially expression analysis UUOVH vs UUOIH/', "/MAplot_cluster_rawcount_", cc, ".pdf", sep = ""))
  plotMA(res, main = paste("cell_", cc, sep = ""))
  dev.off()
}

getwd()
