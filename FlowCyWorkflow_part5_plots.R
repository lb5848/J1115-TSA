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
library(scales)
library(gplots)
library(svglite)

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

outputDirectory <- getwd()
outputDirectory <- paste(outputDirectory, "output", sep = "/")
dir.create(outputDirectory)
setwd(outputDirectory)
load("workspaceFinal.rds")
sce
CATALYST::pbMDS(sce, by = "sample_id", color_by = "condition", features = type_markers(sce), fun = "median", label_by = NULL)

getwd()
plotFolder <- paste(outputDirectory, "CancerDiscoveryPlots_Fig3", sep = "/")
dir.create(plotFolder)
setwd(plotFolder)

# plotAbundances w/ stats

stat.test <- as_tibble(da)
# remove C8 - aggregates and <100 cells/1%
p.adj.signif <- c("ns", rep("***", 3), "ns", rep("***", 3))
group1 <- (rep("baseline", nrow(stat.test)))
group2 <- (rep("more4", nrow(stat.test)))
plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
y.position <- c(90, 3, 55, 9, 50, 2, 15, 25)
stat.test <- cbind(stat.test, group1, group2, p.adj.signif, y.position)
bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
bxp <- bxp + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, size = 2.5)
bxp

ggsave("abundances_stat.svg", plot = last_plot(), dpi = 300)

display.brewer.all(colorblindFriendly = TRUE)

# MDS plot

CATALYST::pbMDS(sce, by = "sample_id", color_by = "condition", 
                features = type_markers(sce), fun = "median", label_by = NULL)
ggsave("MDSplot.svg", plot = last_plot(), dpi = 300)



# UMAP color_by = Clusters facet_by = condition

CATALYST::plotDR(sce, dr = "UMAP", color_by = "cluster_annotation", facet_by = "condition") + 
  scale_color_brewer(palette = "Dark2") +
  geom_density2d(binwidth = 0.001, colour = "black")
ggsave("UMAP_w_contours_bycondition.svg", plot = last_plot(), dpi = 300)

plotAbundances(sce, k = "cluster_annotation", by = "sample_id", group_by = "condition")
ggsave("abubdances_allPts.svg", plot = last_plot(), dpi = 300)

# plotDiffHeatmap
plotDiffHeatmap(sce, da, top_n = 8, all = TRUE, fdr = FDR_cutoff)
ggsave("plotDiffHeatmap.svg", plot = last_plot(), dpi = 300)
