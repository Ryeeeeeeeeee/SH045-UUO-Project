rm(list = ls())
gc()
library(SCENIC)
library(RcisTarget)
library(GENIE3)
library(AUCell)
library(SCopeLoomR)
library(Seurat)
scRNA <- readRDS('F:/kidney/output/UUOVH_vs_UUOIH/UUOVH_vs_UUOIH_merge_DF_SCTharmony_celldenifition_allsubsets.rds')
Idents(scRNA)
scRNAsub <- subset(scRNA,idents = 'Fibroblasts')
exprMat <- as.matrix(scRNAsub@assays$RNA@counts)
cellInfo <- as.data.frame(scRNAsub@meta.data)
data(list="motifAnnotations_mgi_v9", package="RcisTarget")
motifAnnotations_mgi <- motifAnnotations_mgi_v9
scenicOptions <- initializeScenic(org="mgi",
                                  dbDir="D:/scrna/my script/cisTarget_databases/mouse/mm9/", nCores=8)
scenicOptions@inputDatasetInfo$cellInfo <- "cellInfo.Rds"
saveRDS(scenicOptions, file="scenicOptions.Rds")
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept,]
dim(exprMat_filtered)
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
mymethod <- 'R'
library(reticulate)
if(mymethod=='R'){
 
  runGenie3(exprMat_filtered_log, scenicOptions)
}else{
  
  arb.algo = import('arboreto.algo')
  tf_names = getDbTfs(scenicOptions)
  tf_names = Seurat::CaseMatch(
    search = tf_names,
    match = rownames(exprMat_filtered))
  adj = arb.algo$grnboost2(
    as.data.frame(t(as.matrix(exprMat_filtered))),
    tf_names=tf_names, seed=2023L
  )
  colnames(adj) = c('TF','Target','weight')
  saveRDS(adj,file=getIntName(scenicOptions,
                              'genie3ll'))
}
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 3
scenicOptions@settings$seed <- 123
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)

scenicOptions <- readRDS("int/scenicOptions.Rds")
minGenes=20
coexMethod=NULL
vignette("detailedStep_2_createRegulons", package="SCENIC")
nCores <- getSettings(scenicOptions, "nCores")
nCores
tfModules_asDF <- loadInt(scenicOptions, "tfModules_asDF")
head(tfModules_asDF)
if(!is.null(coexMethod)) tfModules_asDF <- tfModules_asDF[which(tfModules_asDF$method %in% coexMethod),]
if(nrow(tfModules_asDF)==0) stop("The co-expression modules are empty.")
if("BiocParallel" %in% installed.packages()) library(BiocParallel); register(MulticoreParam(nCores), default=TRUE) 
msg <- paste0(format(Sys.time(), "%H:%M"), "\tStep 2. Identifying regulons")
if(getSettings(scenicOptions, "verbose")) message(msg)
if(is.na(getDatasetInfo(scenicOptions, "org"))) stop('Please provide an organism (scenicOptions@inputDatasetInfo$org).')

motifAnnot <- getDbAnnotations(scenicOptions)
motifAnnot[1:4,1:4]
if(is.null(names(getSettings(scenicOptions, "dbs")))) 
{
  names(scenicOptions@settings$"dbs") <- scenicOptions@settings$"dbs"
  tmp <- sapply(strsplit(getSettings(scenicOptions, "dbs"),"-", fixed=T), function(x) x[grep("bp|kb",x)])
  if(all(lengths(tmp)>0)) names(scenicOptions@settings$"dbs") <- tmp
}

loadAttempt <- sapply(getDatabases(scenicOptions), dbLoadingAttempt)
if(any(!loadAttempt)) stop("It is not possible to load the following databses: \n",
                           paste(dbs[which(!loadAttempt)], collapse="\n"))

genesInDb <- unique(unlist(lapply(getDatabases(scenicOptions), function(x)
  names(feather::feather_metadata(x)[["types"]]))))
genesInDb[1:10]

tfModules_asDF$TF <- as.character(tfModules_asDF$TF)
tfModules_asDF$Target <- as.character(tfModules_asDF$Target)
allTFs <- getDbTfs(scenicOptions)
allTFs[1:10]
tfModules_asDF <- tfModules_asDF[which(tfModules_asDF$TF %in% allTFs),]
geneInDb <- tfModules_asDF$Target %in% genesInDb
geneInDb[1:10]
missingGene <- sort(unique(tfModules_asDF[which(!geneInDb),"Target"]))
missingGene
if(length(missingGene)>0) 
  warning(paste0("Genes in co-expression modules not available in RcisTargetDatabases: ", 
                 paste(missingGene, collapse=", ")))
tfModules_asDF <- tfModules_asDF[which(geneInDb),]
tfModules_Selected <- tfModules_asDF[which(tfModules_asDF$corr==1),]
# Add a column with the geneSet name (TF_method)
tfModules_Selected <- cbind(tfModules_Selected, geneSetName=paste(tfModules_Selected$TF, tfModules_Selected$method, sep="_"))
head(tfModules_Selected)
tfModules_Selected$geneSetName <- factor(as.character(tfModules_Selected$geneSetName))
allGenes <- unique(tfModules_Selected$Target)
tfModules <- split(tfModules_Selected$Target, tfModules_Selected$geneSetName)

# Add TF to the gene set (used in the following steps, careful if editing)
tfModules <- setNames(lapply(names(tfModules), function(gsn) {
  tf <- strsplit(gsn, "_")[[1]][1]
  unique(c(tf, tfModules[[gsn]]))
}), names(tfModules))
tfModules <- tfModules[which(lengths(tfModules)>=minGenes)]
saveRDS(tfModules, file=getIntName(scenicOptions, "tfModules_forEnrichment"))
print(getIntName(scenicOptions, "tfModules_forEnrichment"))

if(getSettings(scenicOptions, "verbose")) {
  tfModulesSummary <- t(sapply(strsplit(names(tfModules), "_"), function(x) x[1:2]))
  message("tfModulesSummary:")
  print(sort(table(tfModulesSummary[,2])))
}

msg <- paste0(format(Sys.time(), "%H:%M"), "\tRcisTarget: Calculating AUC")
if(getSettings(scenicOptions, "verbose")) message(msg)
motifs_AUC <- lapply(getDatabases(scenicOptions), function(rnkName) {
  ranking <- importRankings(rnkName, columns=allGenes)
  message("Scoring database: ", ranking@description)
  RcisTarget::calcAUC(tfModules, ranking, aucMaxRank=0.03*getNumColsInDB(ranking), nCores=nCores, verbose=FALSE)})
saveRDS(motifs_AUC, file=getIntName(scenicOptions, "motifs_AUC"))

msg <- paste0(format(Sys.time(), "%H:%M"), "\tRcisTarget: Adding motif annotation")
message(msg)
motifEnrichment <- lapply(motifs_AUC, function(aucOutput)
{
  # Extract the TF of the gene-set name (i.e. MITF_w001):
  tf <- sapply(setNames(strsplit(rownames(aucOutput), "_"), rownames(aucOutput)), function(x) x[[1]])
  
  # Calculate NES and add motif annotation (provide tf in 'highlightTFs'):
  addMotifAnnotation(aucOutput, #AUCell score
                     nesThreshold=3, 
                     digits=3, 
                     motifAnnot=motifAnnot,
                     motifAnnot_highConfCat=c("directAnnotation", "inferredBy_Orthology"),
                     motifAnnot_lowConfCat=c("inferredBy_MotifSimilarity",
                                             "inferredBy_MotifSimilarity_n_Orthology"), 
                     highlightTFs=tf
  )
})

motifEnrichment <- do.call(rbind, lapply(names(motifEnrichment), function(dbName){
  cbind(motifDb=dbName, motifEnrichment[[dbName]])
}))
head(motifEnrichment)
saveRDS(motifEnrichment, file=getIntName(scenicOptions, "motifEnrichment_full"))
msg <- paste0("Number of motifs in the initial enrichment: ", nrow(motifEnrichment))
if(getSettings(scenicOptions, "verbose")) message(msg)


motifEnrichment_selfMotifs <- motifEnrichment[which(motifEnrichment$TFinDB != ""),, drop=FALSE]
msg <- paste0("Number of motifs annotated to the corresponding TF: ", nrow(motifEnrichment_selfMotifs))
if(getSettings(scenicOptions, "verbose")) message(msg)
if(nrow(motifEnrichment_selfMotifs)==0) 
  stop("None of the co-expression modules present enrichment of the TF motif: There are no regulons.")

msg <- paste0(format(Sys.time(), "%H:%M"), "\tRcisTarget: Prunning targets")
if(getSettings(scenicOptions, "verbose")) message(msg)
dbNames <- getDatabases(scenicOptions)
motifEnrichment_selfMotifs_wGenes <- lapply(names(dbNames), function(motifDbName){
  ranking <- importRankings(dbNames[motifDbName], columns=allGenes)
  addSignificantGenes(resultsTable=motifEnrichment_selfMotifs[motifEnrichment_selfMotifs$motifDb==motifDbName,],
                      geneSets=tfModules,
                      rankings=ranking,
                      maxRank=5000, method="aprox", nCores=nCores)
})
suppressPackageStartupMessages(library(data.table))
motifEnrichment_selfMotifs_wGenes <- rbindlist(motifEnrichment_selfMotifs_wGenes)
saveRDS(motifEnrichment_selfMotifs_wGenes, file=getIntName(scenicOptions, "motifEnrichment_selfMotifs_wGenes"))

if(getSettings(scenicOptions, "verbose")) 
{
 
  message(format(Sys.time(), "%H:%M"), "\tNumber of motifs that support the regulons: ", nrow(motifEnrichment_selfMotifs_wGenes))
  motifEnrichment_selfMotifs_wGenes[order(motifEnrichment_selfMotifs_wGenes$NES,decreasing=TRUE),][1:5,(1:ncol(motifEnrichment_selfMotifs_wGenes)-1), with=F] 
}

if(!file.exists("output")) dir.create("output")
write.table(motifEnrichment_selfMotifs_wGenes, file=getOutName(scenicOptions, "s2_motifEnrichment"),
            sep="\t", quote=FALSE, row.names=FALSE)


if("DT" %in% installed.packages() && nrow(motifEnrichment_selfMotifs_wGenes)>0)
{
  nvm <- tryCatch({
    colsToShow <- c("motifDb", "logo", "NES", "geneSet", "TF_highConf", "TF_lowConf")
    motifEnrichment_2html <- viewMotifs(motifEnrichment_selfMotifs_wGenes, colsToShow=colsToShow, options=list(pageLength=100))
    
    fileName <- getOutName(scenicOptions, "s2_motifEnrichmentHtml")
    
    dirName <- dirname(fileName)
    fileName <- basename(fileName)
    suppressWarnings(DT::saveWidget(motifEnrichment_2html, fileName))
    file.rename(fileName, file.path(dirName, fileName))
    if(getSettings(scenicOptions, "verbose")) message("Preview of motif enrichment saved as: ", file.path(dirName, fileName))
  }, error = function(e) print(e$message))
}
motifEnrichment.asIncidList <- apply(motifEnrichment_selfMotifs_wGenes, 1, function(oneMotifRow) {
  genes <- strsplit(oneMotifRow["enrichedGenes"], ";")[[1]]
  oneMotifRow <- data.frame(rbind(oneMotifRow), stringsAsFactors=FALSE)
  data.frame(oneMotifRow[rep(1, length(genes)),c("NES", "motif", "highlightedTFs", "TFinDB")], genes, stringsAsFactors = FALSE)
})
class(motifEnrichment.asIncidList)
motifEnrichment.asIncidList <- rbindlist(motifEnrichment.asIncidList)
head(motifEnrichment.asIncidList)
colnames(motifEnrichment.asIncidList) <- c("NES", "motif", "TF", "annot", "gene")
motifEnrichment.asIncidList <- data.frame(motifEnrichment.asIncidList, stringsAsFactors = FALSE)
regulonTargetsInfo <- lapply(split(motifEnrichment.asIncidList, motifEnrichment.asIncidList$TF),
                             function(tfTargets){
                               print(unique(tfTargets$TF))
                               tfTable <- as.data.frame(do.call(rbind, lapply(split(tfTargets, tfTargets$gene), function(enrOneGene){
                                 highConfAnnot <- "**" %in% enrOneGene$annot
                                 enrOneGeneByAnnot <- enrOneGene
                                 if(highConfAnnot) enrOneGeneByAnnot <- enrOneGeneByAnnot[which(enrOneGene$annot == "**"),]
                                 bestMotif <- which.max(enrOneGeneByAnnot$NES)
                                 cbind(TF=unique(enrOneGene$TF), gene=unique(enrOneGene$gene), nMotifs=nrow(enrOneGene),
                                       bestMotif=as.character(enrOneGeneByAnnot[bestMotif,"motif"]), NES=as.numeric(enrOneGeneByAnnot[bestMotif,"NES"]),
                                       highConfAnnot=highConfAnnot)
                               })), stringsAsFactors=FALSE)
                               tfTable[order(tfTable$NES, decreasing = TRUE),]
                             })
regulonTargetsInfo <- rbindlist(regulonTargetsInfo)
colnames(regulonTargetsInfo) <- c("TF", "gene", "nMotifs", "bestMotif", "NES", "highConfAnnot")
head(regulonTargetsInfo)


linkList <- loadInt(scenicOptions, "genie3ll", ifNotExists="null")
if(!is.null(linkList) & ("weight" %in% colnames(linkList)))
{
  if(is.data.table(linkList)) linkList <- as.data.frame(linkList)
  
  uniquePairs <- nrow(unique(linkList[,c("TF", "Target")]))
  if(uniquePairs == nrow(linkList)) {
    linkList <- linkList[which(linkList$weight>=getSettings(scenicOptions, "modules/weightThreshold")),]  
    rownames(linkList) <- paste(linkList$TF, linkList$Target,sep="__")
    regulonTargetsInfo <- cbind(regulonTargetsInfo, Genie3Weight=linkList[paste(regulonTargetsInfo$TF, regulonTargetsInfo$gene,sep="__"),"weight"])
  }else {
    warning("There are duplicated regulator-target (gene id/name) pairs in the co-expression link list.",
            "\nThe co-expression weight was not added to the regulonTargetsInfo table.")
  }
}else warning("It was not possible to add the weight to the regulonTargetsInfo table.")

saveRDS(regulonTargetsInfo, file=getIntName(scenicOptions, "regulonTargetsInfo"))

write.table(regulonTargetsInfo, file=getOutName(scenicOptions, "s2_regulonTargetsInfo"),
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
regulonTargetsInfo_splitByAnnot <- split(regulonTargetsInfo, regulonTargetsInfo$highConfAnnot)
regulons <- NULL
if(!is.null(regulonTargetsInfo_splitByAnnot[["TRUE"]]))
{
  regulons <- lapply(split(regulonTargetsInfo_splitByAnnot[["TRUE"]], regulonTargetsInfo_splitByAnnot[["TRUE"]][,"TF"]), function(x) sort(as.character(unlist(x[,"gene"]))))
}
regulons_extended <- NULL
if(!is.null(regulonTargetsInfo_splitByAnnot[["FALSE"]]))
{
  regulons_extended <- lapply(split(regulonTargetsInfo_splitByAnnot[["FALSE"]],regulonTargetsInfo_splitByAnnot[["FALSE"]][,"TF"]), function(x) unname(unlist(x[,"gene"])))
  regulons_extended <- setNames(lapply(names(regulons_extended), function(tf) sort(unique(c(regulons[[tf]], unlist(regulons_extended[[tf]]))))), names(regulons_extended))
  names(regulons_extended) <- paste(names(regulons_extended), "_extended", sep="")
}
regulons <- c(regulons, regulons_extended)
saveRDS(regulons, file=getIntName(scenicOptions, "regulons"))

incidList <- reshape2::melt(regulons)
head(incidList)
incidMat <- table(incidList[,2], incidList[,1])
saveRDS(incidMat, file=getIntName(scenicOptions, "regulons_incidMat"))

if(getSettings(scenicOptions, "verbose")) 
{
  # Number of regulons and summary of sizes:
  length(regulons) 
  summary(lengths(regulons))
}
scenicOptions@status$current <- 2

scenicOptions@status$current


scenicOptions <- readRDS("int/scenicOptions.rds")
skipBinaryThresholds=FALSE # Whether to skip the automatic binarization step
skipHeatmap=FALSE # hether to plot the AUC heatmap
skipTsne=FALSE # Whether to plot the t-SNE

regulons <- loadInt(scenicOptions, "regulons")
regulons <- regulons[order(lengths(regulons), decreasing=TRUE)]
regulons <- regulons[lengths(regulons)>=10]
if(length(regulons) <2)  stop("Not enough regulons with at least 10 genes.")
# Add the TF to the regulon (keeping it only once) & rename regulon
regulons <- setNames(lapply(names(regulons), function(tf) sort(unique(c(gsub("_extended", "", tf), regulons[[tf]])))), names(regulons))
names(regulons) <- paste(names(regulons), " (",lengths(regulons), "g)", sep="")
saveRDS(regulons, file=getIntName(scenicOptions, "aucell_regulons"))
msg <- paste0(format(Sys.time(), "%H:%M"), "\tStep 3. Analyzing the network activity in each individual cell")
if(getSettings(scenicOptions, "verbose")) message(msg)

msg <- paste0("\nNumber of regulons to evaluate on cells: ", length(regulons),
              "\nBiggest (non-extended) regulons: \n",
              paste("\t", grep("_extended",names(regulons),invert = T, value = T)[1:10], collapse="\n")) 
if(getSettings(scenicOptions, "verbose")) message(msg)
#AUCell

if(is.data.frame(exprMat)) 
{
  supportedClasses <- paste(gsub("AUCell_buildRankings,", "", methods("AUCell_buildRankings")), collapse=", ")
  supportedClasses <- gsub("-method", "", supportedClasses)
  
  stop("'exprMat' should be one of the following classes: ", supportedClasses, 
       "\n(data.frames are not supported. Please, convert the expression matrix to one of these classes.)")
}

set.seed(getSettings(scenicOptions,"seed"))
tryCatch({
  pdf(file=paste0(getIntName(scenicOptions, "aucell_genesStatsPlot"),'.my.pdf'))
           
  aucellRankings <- AUCell_buildRankings(exprMat, nCores=nCores, 
                                         plotStats=TRUE, verbose=getSettings(scenicOptions, "verbose"))
  abline(v=aucellRankings@nGenesDetected["1%"], col="skyblue3", lwd=5, lty=3)
  dev.off()
},error = function(e) {
  message("Catched error in AUCell_buildRankings() or in the histogram plot: ", e$message)
})
saveRDS(aucellRankings, file=getIntName(scenicOptions, "aucell_rankings"))

regulonAUC <- AUCell_calcAUC(regulons, aucellRankings, 
                             aucMaxRank=aucellRankings@nGenesDetected["1%"], nCores=nCores)

# Order the modules by similarity, for easier exploration in the upcoming steps & save
regulonOrder <- orderAUC(regulonAUC) 
regulonAUC <- regulonAUC[regulonOrder,]
saveRDS(regulonAUC, file=getIntName(scenicOptions, "aucell_regulonAUC"))

cells_AUCellThresholds <- NULL
if(!skipBinaryThresholds)
{
  cells_AUCellThresholds <- AUCell_exploreThresholds(regulonAUC, 
                                                     smallestPopPercent=getSettings(scenicOptions,"aucell/smallestPopPercent"),
                                                     assignCells=TRUE, plotHist=FALSE, 
                                                     verbose=FALSE, nCores=nCores)
  saveRDS(cells_AUCellThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
  
  # Get cells assigned to each regulon
  regulonsCells <- getAssignments(cells_AUCellThresholds)
  
  ### Save threshold info as text (e.g. to edit/modify...)
  trhAssignment <- getThresholdSelected(cells_AUCellThresholds)
  trhAssignment <- signif(trhAssignment, 3) 
  commentsThresholds <- sapply(cells_AUCellThresholds, function(x) unname(x$aucThr$comment))
  
  table2edit <- cbind(regulon=names(cells_AUCellThresholds),
                      threshold=trhAssignment[names(cells_AUCellThresholds)],
                      nCellsAssigned=lengths(regulonsCells)[names(cells_AUCellThresholds)],
                      AUCellComment=commentsThresholds[names(cells_AUCellThresholds)],
                      nGenes=gsub("[\\(g\\)]", "", regmatches(names(cells_AUCellThresholds), gregexpr("\\(.*?\\)", names(cells_AUCellThresholds)))),
                      clusteringOrder=1:length(cells_AUCellThresholds),
                      clusterGroup=regulonClusters[names(cells_AUCellThresholds)],
                      onlyNonDuplicatedExtended=(names(cells_AUCellThresholds) %in% onlyNonDuplicatedExtended(names(cells_AUCellThresholds))),
                      personalNotes="")
  write.table(table2edit, file=getIntName(scenicOptions, "aucell_thresholdsTxt"), row.names=F, quote=F, sep="\t")
  rm(trhAssignment)
}

if(!skipHeatmap){
  nCellsHeatmap <- min(500, ncol(regulonAUC))
  cells2plot <- sample(colnames(regulonAUC), nCellsHeatmap)
  
  cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null")  
  if(!is.null(cellInfo)) cellInfo <- data.frame(cellInfo)[cells2plot,,drop=F]
  colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars"), ifNotExists="null")
  
  fileName <- getOutName(scenicOptions, "s3_AUCheatmap")
  
  fileName <- .openDevHeatmap(fileName=fileName, devType=getSettings(scenicOptions, "devType"))
  NMF::aheatmap(getAUC(regulonAUC)[,cells2plot],
                annCol=cellInfo,
                annColor=colVars,
                main="AUC",
                sub=paste("Subset of",nCellsHeatmap," random cells"),
                filename=paste0(fileName,'.pdf'))
}
#tsne
if(!skipTsne){
  tSNE_fileName <- tsneAUC(scenicOptions, aucType="AUC", onlyHighConf=FALSE) # default: nPcs, perpl, seed, tsne prefix
  tSNE <- readRDS(tSNE_fileName)
  
  # AUCell (activity) plots with the default tsne, as html: 
  fileName <- getOutName(scenicOptions, "s3_AUCtSNE_colAct")
  plotTsne_regulonActivityHTML(scenicOptions, exprMat, fileName, tSNE) #open the resulting html locally
  
  # Plot cell properties:
  sub <- ""; if("type" %in% names(tSNE)) sub <- paste0("t-SNE on ", tSNE$type)
  cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null") 
  colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars"), ifNotExists="null")
  pdf(paste0(getOutName(scenicOptions, "s3_AUCtSNE_colProps"),".pdf"))
  plotTsne_cellProps(tSNE$Y, cellInfo=cellInfo, colVars=colVars, cex=1, sub=sub)
  dev.off()
}
# Finished. Update status.
object@status$current <- 3
saveRDS(scenicOptions,file = 'scenicOptions.rds')


exprMat_log <- log2(exprMat+1)
aucellApp <- plotTsne_AUCellApp(scenicOptions,exprMat_log)
######
skipBoxplot = FALSE
skipHeatmaps = FALSE
skipTsne = FALSE
exprMat = NULL
nCores <- getSettings(scenicOptions, "nCores")
regulonAUC <- tryCatch(loadInt(scenicOptions, "aucell_regulonAUC"), 
                       error = function(e) {
                         if (getStatus(scenicOptions, asID = TRUE) < 3) 
                           e$message <- paste0("It seems the regulons have not been scored on the cells yet. Please, run runSCENIC_3_scoreCells() first.\n", 
                                               e$message)
                         stop(e)
                       })
thresholds <- loadInt(scenicOptions, "aucell_thresholds")
thresholds <- getThresholdSelected(thresholds)
regulonsCells <- setNames(lapply(names(thresholds), function(x) {
  trh <- thresholds[x]
  names(which(getAUC(regulonAUC)[x, ] > trh))
}), names(thresholds))
regulonActivity <- reshape2::melt(regulonsCells)
binaryRegulonActivity <- t(table(regulonActivity[, 1], regulonActivity[, 
                                                                       2]))
class(binaryRegulonActivity) <- "matrix"
saveRDS(binaryRegulonActivity, file = getIntName(scenicOptions, 
                                                 "aucell_binary_full"))
if (nrow(binaryRegulonActivity) == 0) 
  stop("No cells passed the binarization.")
binaryRegulonActivity_nonDupl <- binaryRegulonActivity[which(rownames(binaryRegulonActivity) %in% 
                                                               onlyNonDuplicatedExtended(rownames(binaryRegulonActivity))), 
]
saveRDS(binaryRegulonActivity_nonDupl, file = getIntName(scenicOptions, 
                                                         "aucell_binary_nonDupl"))
minCells <- ncol(binaryRegulonActivity) * 0.01
msg <- paste0("Binary regulon activity: ", nrow(binaryRegulonActivity_nonDupl), 
              " TF regulons x ", ncol(binaryRegulonActivity), " cells.\n(", 
              nrow(binaryRegulonActivity), " regulons including 'extended' versions)\n", 
              sum(rowSums(binaryRegulonActivity_nonDupl) > minCells), 
              " regulons are active in more than 1% (", minCells, 
              ") cells.")
if (getSettings(scenicOptions, "verbose")) 
  message(msg)
if (!skipBoxplot) {
  .openDev(fileName = getOutName(scenicOptions, "s4_boxplotBinaryActivity"), 
           devType = getSettings(scenicOptions, "devType"))
  par(mfrow = c(1, 2))
  boxplot(rowSums(binaryRegulonActivity_nonDupl), main = "nCells per regulon", 
          sub = "number of cells \nthat have the regulon active", 
          col = "darkolivegreen1", border = "#001100", lwd = 2, 
          frame = FALSE)
  boxplot(colSums(binaryRegulonActivity_nonDupl), main = "nRegulons per Cell", 
          sub = "number of regulons \nactive per cell", col = "darkolivegreen1", 
          border = "#001100", lwd = 2, frame = FALSE)
  dev.off()
}
if (!skipHeatmaps) {
  regulonSelection <- loadInt(scenicOptions, "aucell_regulonSelection", 
                              ifNotExists = "null", verbose = FALSE)
  if (is.null(regulonSelection)) 
    regulonSelection <- regulonSelections(binaryRegulonActivity, 
                                          binaryRegulonActivity_nonDupl, minCells)
  cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, 
                                                     "cellInfo"), ifNotExists = "null")
  cellInfo <- data.frame(cellInfo)
  colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, 
                                                    "colVars"), ifNotExists = "null")
  for (selRegs in names(regulonSelection$labels)) {
    if (length(regulonSelection[[selRegs]]) > 0) {
      regulonSelection[[selRegs]] <- regulonSelection[[selRegs]][which(regulonSelection[[selRegs]] %in% 
                                                                         rownames(binaryRegulonActivity))]
      binaryMat <- binaryRegulonActivity[regulonSelection[[selRegs]], 
                                         , drop = FALSE]
      if (nrow(binaryMat) > 0) {
        fileName <- paste0(getOutName(scenicOptions, 
                                      "s4_binaryActivityHeatmap"), selRegs)
        fileName <- .openDevHeatmap(fileName = fileName, 
                                    devType = getSettings(scenicOptions, "devType"))
        rowv <- ifelse(nrow(binaryMat) >= 2, T, NA)
        colv <- ifelse(ncol(binaryMat) >= 2, T, NA)
        NMF::aheatmap(binaryMat, scale = "none", revC = TRUE, 
                      main = selRegs, annCol = cellInfo[colnames(binaryMat), 
                                                        , drop = FALSE], annColor = colVars, Rowv = rowv, 
                      Colv = colv, color = c("white", "black"), 
                      filename = fileName)
        if (getSettings(scenicOptions, "devType") != 
            "pdf") 
          dev.off()
      }
      else {
        if (getSettings(scenicOptions, "verbose")) 
          message(paste0("No regulons to plot for regulon selection '", 
                         selRegs, "'. Skipping."))
      }
    }
  }
}
if (!skipTsne) {
  tSNE_fileName <- tsneAUC(scenicOptions, aucType = "Binary", 
                           filePrefix = getIntName(scenicOptions, "tsne_prefix"), 
                           onlyHighConf = FALSE)
  if (!is.null(tSNE_fileName)) {
    tSNE <- readRDS(tSNE_fileName)
    fileName <- getOutName(scenicOptions, "s4_binarytSNE_colAct")
    plotTsne_AUCellHtml(scenicOptions, exprMat, fileName, 
                        tSNE)
    sub <- ""
    if ("type" %in% names(tSNE)) 
      sub <- paste0("t-SNE on ", tSNE$type)
    cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, 
                                                       "cellInfo"), ifNotExists = "null")
    colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, 
                                                      "colVars"), ifNotExists = "null")
    pdf(paste0(getOutName(scenicOptions, "s4_binarytSNE_colProps"), 
               ".pdf"))
    plotTsne_cellProps(tSNE$Y, cellInfo = cellInfo, 
                       colVars = colVars, cex = 1, sub = sub)
    dev.off()
  }
}
scenicOptions@status$current <- 4
invisible(scenicOptions)

nPcs <- c(5,15,50)
nperpl <- c(5,15,50)
scenicOptions@settings$seed <- 123
fileNames <- tsneAUC(scenicOptions,aucType = "AUC",
                     nPcs = nPcs,
                     perpl = nperpl)

fileNames
fileNames2 <- tsneAUC(scenicOptions,aucType = "AUC",
                      nPcs = nPcs,perpl = nperpl,
                      onlyHighConf = TRUE,
                      filePrefix = 'int/tSNE_oHC')
fileNames2
fileNames3 <- tsneAUC(scenicOptions,aucType = "Binary",
                      nPcs = nPcs,perpl = nperpl,
                      onlyHighConf = TRUE,
                      filePrefix = 'int/tSNE_oHC')
fileNames3
fileNames <- paste0("int/",
                    grep(".Rds",
                         grep("tSNE_AUC", list.files("int"),
                              value = T,perl = T),value = T))

plotTsne_compareSettings(fileNames, scenicOptions, 
                         showLegend=FALSE, varName = "CellType",
                         cex=.5)
par(mfrow=c(3,3))
fileNames <- paste0("int/",
                    grep(".Rds",
                         grep("tSNE_AUC", list.files("int"),
                              value = T,perl = T),value = T))
plotTsne_compareSettings(fileNames, scenicOptions, 
                         showLegend=FALSE, varName = "CellType",
                         cex=.5)


scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 5
scenicOptions@settings$defaultTsne$perpl <- 15
saveRDS(scenicOptions,file = "scenicOptions.Rds")


print(tsneFileName(scenicOptions))
tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions,'aucell_regulonAUC')

par(mfrow=c(2,3))

AUCell::AUCell_plotTSNE(
  tSNE_scenic$Y,
  exprMat_log,
     aucell_regulonAUC[
       onlyNonDuplicatedExtended(
         rownames(aucell_regulonAUC))[
           c('Dlx5','Sox10','Sox9')
         ],
     ],
  plots = 'Expression')

Cairo::CairoPDF('output/Step4_con_RegulonActivity_tSNE_colByAUC.pdf',
                width=20,height = 15)
par(mfrow=c(4,6))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y,cellsAUC = aucell_regulonAUC,plots = "AUC")
dev.off()

par(mfrow=c(1,2))
regulonMames <- c('Dlx5','Sox10')
SCENIC::plotEmb_rgb(scenicOptions,
                    regulonNames = regulonMames,aucType = 'AUC',
                      aucMaxContrast = 0.6,offColor = 'lightgray')
regulonMames <- list(red=c('Sox10','Sox8'),
                     green=c('Irf1'),
                     blue=C('Tef'))
SCENIC::plotEmb_rgb(scenicOptions,
                    regulonNames = regulonMames,aucType = 'Binary',
                    aucMaxContrast = 0.6,offColor = 'lightgray')

regulons <- loadInt(scenicOptions,'regulons')
regulons
regulonTargetsInfo<- loadInt(scenicOptions,'regulonTargetsInfo')
tableSubset <- regulonTargetsInfo[TF=='Stat6' & highConfAnnot==TRUE]
viewMotifs(tableSubset,options = list(pageLength=5))

motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions,'motifEnrichment_selfMotifs_wGenes')
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=='Dlx5']
viewMotifs(tableSubset)

cellInfo <- data.frame(seuratCluster=Idents(seuratObject))
loomPath <- system.file(package="SCENIC", "examples/mouseBrain_toy.loom")

loom <- open_loom(loomPath)
cellInfo <- get_cell_annotation(loom)
regulonAUC <- loadInt(scenicOptions,'aucell_regulonAUC')
cellInfo <- cellInfo[colnames(regulonAUC),]
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo),cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,as.character(cells)]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType),center=T,scale = T))
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name = 'Regulon activity')

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c('Regulon','CellType','RelativeActivity')
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)

minPerc <- 0.7
binaryRegulonActivity <- loadInt(scenicOptions,'aucell_binary_nonDipl')
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in%colnames(binaryRegulonActivity)),,drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells),cellInfo_binarizedCells$CellType),
                                               function(cells) rowMeans(binaryRegulonActivity[,cells,drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
ComplexHeatmap::Heatmap(binaryActPerc_subset,
                        name='Regulon activity (%)', col =c('white','pink','red'))
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Binarized,
                        name='Regulon activity (%)', col =c('white','pink','red'))

rss <- calcRSS(AUC = getAUC(regulonAUC),
               cellAnnotation = cellInfo[colnames(regulonAUC),'CellType'])
rssPlot <- plotRSS(rss)

plotly::ggplotly(rssPlot$plot)
plotRSS_oneSet(rss,setName = 'interneurons')
