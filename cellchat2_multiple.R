Here is the provided code with added English comments for each step:
  
  ```r
# Load necessary libraries
library(CellChat)
library(patchwork)
library(cowplot)

# Clear the environment
rm(list = ls()) 
gc()

# Load the pre-saved CellChat objects
cellchat.UUOVH <- readRDS("F:/kidney/output/UUOVH/cellchat_forinfl/UUOVH_cellchat.rds")
cellchat.UUOIH <- readRDS("F:/kidney/output/UUOIH/cellchat_forinfl/UUOIH_cellchat_removeF2.rds")

# Combine CellChat objects into a list
object.list <- list(UUOVH = cellchat.UUOVH, UUOIH = cellchat.UUOIH)

# Merge the CellChat objects
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

# Compare interactions between groups (1 and 2)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

# Rank and visualize network interactions
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

# Visualize differences in interactions
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

# Visualize circle plots for each group in the object list
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

# Visualize heatmaps
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

# Visualize circle plots for signaling pathways
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

# Visualize heatmaps for signaling pathways
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, 
                               color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

# Define cell types and merge interactions
levels(object.list[[1]]@idents) 
levels(object.list[[2]]@idents) 
group.cellType <- c("Kidney_cells","Endo_CapillaryArterial","Endo_LargeArtery","Endo_Unknown","Endo_Vein","F0","F1","F3","F4","F5","Kidney_cells","Macrophages","Kidney_cells","Dendritic_cells","T_cells","Unknown_infs","B_cells","Neutrophils","Kidney_cells","Myeloid_cells","Kidney_cells","Kidney_cells")
group.cellType <- factor(group.cellType, levels = c("Kidney_cells","Endo_CapillaryArterial","Endo_LargeArtery","Endo_Unknown","Endo_Vein","F0","F1","F3","F4","F5","Macrophages","Dendritic_cells","T_cells","Dendritic_cells","B_cells","Neutrophils","Myeloid_cells"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})

# Merge CellChat objects again after defining cell types
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Visualize differences in interactions after merging
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = F)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = F)

# Visualize circle plots for merged interactions
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= F, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

# Visualize signaling roles
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) +
    colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]],
                                               title = names(object.list)[i], 
                                               weight.MinMax = weight.MinMax)
}

# Customize and save signaling role scatter plots
gg[[1]] <- netAnalysis_signalingRole_scatter(object.list[[1]], title = names(object.list)[1], weight.MinMax = weight.MinMax, dot.alpha=50, dot.size = 0.5)+ scale_y_continuous(limits = c(0,7.5))+ scale_x_continuous(limits = c(0,15))+ geom_point(shape = 18, size =3 )
gg[[2]] <- netAnalysis_signalingRole_scatter(object.list[[2]], title = names(object.list)[2], weight.MinMax = weight.MinMax, dot.alpha=50, dot.size = 0.5)+ scale_y_continuous(limits = c(0,7.5))+ scale_x_continuous(limits = c(0,15))+ geom_point(shape = 18, size =3 )
pdf("E:/output/UUOVH_vs_UUOIH/cellchat_sub_removeF2/single_test.pdf")
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) +
    colSums(x@net$count)-diag(x@net$count)})
num.link

weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
weight.MinMax

# Compute centrality for each object and visualize signaling roles
gg <- list()
for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax, dot.alpha=50, dot.size = 0.5)+ scale_y_continuous(limits = c(0,7.5))+ scale_x_continuous(limits = c(0,15))+ geom_point(shape = 18, size =3 )
}
patchwork::wrap_plots(plots = gg)
dev.off()

# Visualize differential signaling roles
pdf("E:/output/UUOVH_vs_UUOIH/cellchat_sub_removeF2/single_test2.pdf")
netAnalysis_diff_signalingRole_scatter(cellchat, comparison = c(1,2), title = "test")
dev.off()

# Load ComplexHeatmap library
suppressMessages(library(ComplexHeatmap))

# Visualize signaling role heatmaps
i = 1
# Outgoing signaling heatmap
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 15)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

# Incoming signaling heatmap
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 15, color.heatmap = "RdBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway


                    