rm(list = ls())
# Load packages
library(ggrastr)
library(rstudioapi)
library(devtools)
library("flowCore")
library("flowWorkspace")
library(cytofCore)
library(FlowSOM)
# library(cytofkit)
library(cluster)
library(Rtsne)
library(ggplot2)
library(dplyr)
library(flowViz)
library(scales)
library(ggthemes)
library(RColorBrewer)
library(uwot)
library(CATALYST)
library(Rphenograph)
library(diffcyt)
library(SummarizedExperiment)
library(stringr)
library(ggcyto)
library(SingleCellExperiment)
library(scran)
library(scater)
library(flowVS)
library(readxl)
library(flowStats)
library(FlowSOMworkshop)
library(flowAI)
library(CytoNorm)
library(PeacoQC)
library(CytoML)
options(java.parameters="-Xmx60G")
library(tidyverse)
library(data.table)
# library(scuttle)
# library(iMUBAC)
library(ggpubr)




# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory
# Define workingDirectory
wdName <- "200906_Working_DirectoryFCS"
workingDirectory <- paste(PrimaryDirectory, wdName, sep = "/")

setwd(workingDirectory)

# Load workspace and SCEobject

sce <- readRDS("SCE_part2_DR.rds")
CATALYST::pbMDS(sce, by = "sample_id", color_by = "condition", features = type_markers(sce), fun = "median", label_by = NULL)

# choose clustering based on delta_area(sce)
delta_area(sce)
# check abundances
plotAbundances(sce, k = "meta8", by = "cluster_id", group_by = "condition")
plotAbundances(sce, k = "meta10", by = "cluster_id", group_by = "condition")
plotAbundances(sce, k = "meta12", by = "cluster_id", group_by = "condition")
# check ExprHeatmap
plotExprHeatmap(sce, features = type_markers(sce), k = "meta8", by = "cluster_id",
                scale = "last", bars = TRUE, perc = TRUE, fun = "mean")
plotClusterExprs(sce, features = type_markers(sce), k = "meta8")

# cluster merging
CD8pos <- c(4, 8, 6, 7)
CD4pos <- c(1, 5, 2, 3)
DN <- c(4, 1)
DP <- c(7, 3)
IFNgpos <- c(8, 5)
TNFapos <- c(6, 2)
cluster_anno <- c(1:8)
cluster_anno[CD8pos] <- "CD8+"
cluster_anno[CD4pos] <- "CD4+"
cluster_anno[DN] <- paste(cluster_anno[DN], "TNFa- IFNg-", sep = " ")
cluster_anno[DP] <- paste(cluster_anno[DP], "TNFa+ IFNg+", sep = " ")
cluster_anno[IFNgpos] <- paste(cluster_anno[IFNgpos], "TNFa- IFNg+", sep = " ")
cluster_anno[TNFapos] <- paste(cluster_anno[TNFapos], "TNFa+ IFNg-", sep = " ")
merging_table <- cbind(c(1:8), cluster_anno)
merging_table <- as.data.frame(merging_table)
colnames(merging_table)[1] <- "meta8"
merging_table$cluster_anno <- factor(merging_table$cluster_anno)
sce <- mergeClusters(sce, k = "meta8", table = merging_table, id = "cluster_annotation")
plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
# Stat analysis
FDR_cutoff <- 0.05
ei <- sce@metadata$experiment_info


# DA using edgeR
design <- createDesignMatrix(ei,
                             cols_design = c("condition", "patient_id"))

contrast <- createContrast(c(0, 0.25, 0.25, 0.25, 0.25, rep(0, 14)))

nrow(contrast) == ncol(design)

# min_cells = 100, min_samples = 3
out_DA <- diffcyt(sce,
                  experiment_info = ei, design = design, contrast = contrast,
                  analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE,
                  transform = FALSE, normalize = TRUE, min_cells = 100, min_samples = 3)

da <- rowData(out_DA$res)
plotDiffHeatmap(sce, da, top_n = 10, all = TRUE, fdr = FDR_cutoff) 


# Save current workspace
save(list = ls(), file = "workspaceSCEDA.rds")
