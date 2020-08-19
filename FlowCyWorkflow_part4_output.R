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
load("workspaceSCEDA.rds")
sce
CATALYST::pbMDS(sce, by = "sample_id", color_by = "condition", features = type_markers(sce), fun = "median", label_by = NULL)

# Save .fcs per cluster
outputDirectory <- getwd()
outputDirectory <- paste(outputDirectory, "output", sep = "/")
dir.create(outputDirectory)
setwd(outputDirectory)

plotClusterExprs(sce, k = "meta8", features = "type")

# TNFa+ 7,8,6,2
# IFNg+ 7,6
# CD8+ 7,8,4,5
# CD4+ 3,6,1,2
# GrB+ 7,5,3,6
annotation_table <- as.data.frame(cbind(c(1:8), 
                                        c("CD4+ Cytokine-", "CD4+ TNFa+", "CD4+ GrB+", "CD8+ Cytokine-", "CD8+ GrB+", 
                                          "CD4+ TNFa+ IFNg+ GrB+", "CD8+ TNFa+ IFNg+ GrB+", "CD8+ TNFa+")))
colnames(annotation_table) <- c("meta8", "FinalClusters")

annotation_table$FinalClusters <- factor(annotation_table$FinalClusters, 
                                         levels = c("CD4+ Cytokine-", "CD4+ TNFa+", "CD4+ GrB+", "CD4+ TNFa+ IFNg+ GrB+","CD8+ Cytokine-",
                                                    "CD8+ TNFa+", "CD8+ GrB+", "CD8+ TNFa+ IFNg+ GrB+"))
sce <- mergeClusters(sce, k = "meta8", 
                     table = annotation_table, id = "cluster_annotation", overwrite = TRUE)

# keep_dr = TRUE not all cells have DR
sce$cluster_annotation <- cluster_ids(sce, "cluster_annotation")
flowSet <- sce2fcs(sce, split_by = "cluster_annotation", keep_cd = TRUE, keep_dr = FALSE, assay = "counts")
write.flowSet(flowSet, outdir = outputDirectory, filename = "bycluster")
merged <- sce2fcs(sce, split_by = NULL, keep_cd = TRUE, keep_dr = FALSE, assay = "counts")
write.FCS(merged, filename = "merged.fcs")
by_condition <- sce2fcs(sce, split_by = "condition", keep_cd = TRUE, keep_dr = FALSE, assay = "counts")
write.flowSet(by_condition, outdir = outputDirectory, filename = "bycondition")
# Save final workspace
save(list = ls(), file = "workspaceFinal.rds")

