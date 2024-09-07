```r
# Load required libraries
library(Seurat)
library(lsei)

# Load the Seurat objects
adata <- readRDS(file = 'F:/kidney/output/UUOVH_vs_UUOIH/Endothelial/test/Endothelial_Seurat.rds')
adata2 <- readRDS(file = 'F:/kidney/output/kuppe/Endo.rds')

# Set the active identity class to 'Seurat_cluster' for both Seurat objects
Idents(adata) <- 'Seurat_cluster'
Idents(adata2) <- 'Seurat_cluster'

# Find shared genes between the two datasets
shared_gene <- intersect(rownames(adata), rownames(adata2))

# Aggregate expression for the first Seurat object
seu.obj.mat <- Seurat::AggregateExpression(adata, assays = "RNA", features = shared_gene, slot = "counts")$RNA
seu.obj.mat <- seu.obj.mat / rep(colSums(seu.obj.mat), each = nrow(seu.obj.mat))
seu.obj.mat <- log10(seu.obj.mat * 100000 + 1)

# Aggregate expression for the second Seurat object
seu.obj.mat2 <- Seurat::AggregateExpression(adata2, assays = "RNA", features = shared_gene, slot = "counts")$RNA
seu.obj.mat2 <- seu.obj.mat2 / rep(colSums(seu.obj.mat2), each = nrow(seu.obj.mat2))
seu.obj.mat2 <- log10(seu.obj.mat2 * 100000 + 1)

# Initialize lists to store correlation values
list1 <- list()
for (cluster in seq(1, ncol(seu.obj.mat))) {
  
  # Predict a with b
  seu.obj.gene <- seu.obj.mat[, cluster]
  
  gene_fc <- seu.obj.gene / apply(seu.obj.mat, 1, median)
  gene_list1 <- names(sort(gene_fc, decreasing = TRUE)[1:200])
  
  gene_fc <- seu.obj.gene / apply(seu.obj.mat[,-cluster], 1, max)
  gene_list2 <- names(sort(gene_fc, decreasing = TRUE)[1:200])
  
  gene_list <- unique(c(gene_list1, gene_list2))
  
  Ta <- seu.obj.mat[gene_list, cluster]
  Mb <- seu.obj.mat2[gene_list,]
  solv <- nnls(Mb, Ta)
  corr <- solv$x
  
  list1[[cluster]] <- corr
}

list2 <- list()
for (cluster in seq(1, ncol(seu.obj.mat2))) {
  
  # Predict a with b
  seu.obj.gene <- seu.obj.mat2[, cluster]
  
  gene_fc <- seu.obj.gene / apply(seu.obj.mat2, 1, median)
  gene_list1 <- names(sort(gene_fc, decreasing = TRUE)[1:200])
  
  gene_fc <- seu.obj.gene / apply(seu.obj.mat2[,-cluster], 1, max)
  gene_list2 <- names(sort(gene_fc, decreasing = TRUE)[1:200])
  
  gene_list <- unique(c(gene_list1, gene_list2))
  
  Ta <- seu.obj.mat2[gene_list, cluster]
  Mb <- seu.obj.mat[gene_list,]
  solv <- nnls(Mb, Ta)
  corr <- solv$x
  
  list2[[cluster]] <- corr
}

# Combine the results into matrices
mat1 <- do.call(rbind, list1)
mat2 <- do.call(rbind, list2)

# Calculate beta values
beta <- 2 * (mat1 + 0.01) *  t(mat2 + 0.01)

# Set row and column names for the beta matrix
row.names(beta) <- paste0("C", 1:nrow(beta) - 1)
colnames(beta) <- paste0("C", 1:ncol(beta) - 1)

# Plot heatmap of beta values
pheatmap::pheatmap(beta, cluster_rows = FALSE, cluster_cols = FALSE)



























