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
wdName <- "200819_Working_DirectoryFCS"
workingDirectory <- paste(PrimaryDirectory, wdName, sep = "/")

setwd(workingDirectory)

# Load workspace and SCEobject

sce <- readRDS("SCE_part2_DR.rds")
CATALYST::pbMDS(sce, by = "sample_id", color_by = "condition", features = type_markers(sce), fun = "median", label_by = NULL)

# choose clustering based on delta_area(sce)
delta_area(sce)
# check abundances
plotAbundances(sce, k = "meta8", by = "cluster_id", group_by = "condition")
# check ExprHeatmap
plotExprHeatmap(sce, features = type_markers(sce), k = "meta8", by = "cluster_id",
                scale = "first", bars = TRUE, perc = TRUE, fun = "mean")
plotClusterExprs(sce, features = type_markers(sce), k = "meta8")

# Stat analysis
FDR_cutoff <- 0.05
ei <- sce@metadata$experiment_info


# DA using edgeR
design <- createDesignMatrix(ei,
                             cols_design = c("condition", "patient_id"))

contrast <- createContrast(c(0, 1, 1, 1, 1, rep(0, 14)))

nrow(contrast) == ncol(design)

# min_cells = 100, min_samples = 3
out_DA <- diffcyt(sce,
                  experiment_info = ei, design = design, contrast = contrast,
                  analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "meta8", verbose = TRUE, subsampling = TRUE,
                  transform = FALSE, normalize = TRUE, min_cells = 100, min_samples = 3)

da <- rowData(out_DA$res)
plotDiffHeatmap(sce, da, top_n = 10, all = TRUE, fdr = FDR_cutoff) 


# Save current workspace
save(list = ls(), file = "workspaceSCEDA.rds")
