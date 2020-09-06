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
# options(java.parameters="-Xmx60G")
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

# Define flow_files directory
dirName <- "flow_files"
dirDirectory <- paste(PrimaryDirectory, dirName, sep = "/")
dir.create(dirDirectory)


# Define fcs directory
fcsName <- "fcs_files"
fcsDirectory <- paste(dirDirectory, fcsName, sep = "/")
dir.create(fcsDirectory)

# Define csv directory
csvName <- "csv_files"
csvDirectory <- paste(dirDirectory, csvName, sep = "/")
dir.create(csvDirectory)


# Define fcs from csv directory
csv2fcsName <- "fcs_fromCSV_files"
csv2fcsDirectory <- paste(dirDirectory, csv2fcsName, sep = "/")
dir.create(csv2fcsDirectory)

# Define workingDirectory
wdName <- "200906_Working_DirectoryFCS"
workingDirectory <- paste(PrimaryDirectory, wdName, sep = "/")
dir.create(workingDirectory)

# List files

CSVfiles <- list.files(csvDirectory, pattern = ".csv$", full = FALSE)
fileName <- gsub(".csv", ".fcs", CSVfiles)

convertCSVtoFCS <- function(CSVfiles, csvDirectory, csv2fcsDirectory, fileName){
  for(i in c(1:length(CSVfiles))){
    data <- read.csv(paste(csvDirectory, CSVfiles[i], sep = "/"))
    print(CSVfiles[i])
    print(fileName[i])
    cytofCore.write.FCS(as.matrix(data), 
                        filename = paste(csv2fcsDirectory, fileName[i], sep = "/"),
                        what = "numeric")
  }
}

convertCSVtoFCS(CSVfiles = CSVfiles, csvDirectory = csvDirectory, csv2fcsDirectory = csv2fcsDirectory, fileName = fileName)

# read flowSet
# Create flowSet from FCSfiles
FCSfiles <- list.files(csv2fcsDirectory, pattern = ".fcs$", full = FALSE)
flowSet <- read.flowSet(files = FCSfiles, path = csv2fcsDirectory, truncate_max_range = FALSE)
colnames(flowSet)

# Load sample_metadata
sample_mdName <- "sample_metadata.xlsx"
sample_md <- read_excel(path = paste(PrimaryDirectory, sample_mdName, sep = "/"), col_names = TRUE) 
sample_md <- as.data.table(sample_md)
sample_md$patient_id <- paste0("Pt_", sample_md$patient_id)
sample_md$sample_id <- paste(sample_md$condition, sample_md$patient_id, sep = "_")
sample_md
sample_md[, sample_id := factor(sample_id)]
lv <- levels(sample_md$sample_id)
lv <- c(lv[12:21], lv[22:33], lv[1:11], lv[34:41])
sample_md[, sample_id := factor(sample_id, levels = lv)]
sample_md[, condition := factor(condition, levels = c("baseline", "C3D14", "1year", "less3", "more4"))]
colnames(sample_md)

panel_mdName <- "panel_md.xlsx"
panel_md <- read_excel(path = paste(PrimaryDirectory, panel_mdName, sep = "/"), col_names = TRUE)
panel_md <- as.data.table(panel_md)

colname_fs <- panel_md$fcs_colname
length(colnames(flowSet)) == length(colname_fs)
colnames(flowSet) <- colname_fs
colnames(flowSet)

# Preprocessing QC

setwd(workingDirectory)

# QC_dir <- "QC"
# if(!dir.exists(QC_dir)){
#   dir.create(QC_dir)
#   dir.create(file.path(QC_dir, "flowAI"))
#   dir.create(file.path(QC_dir, "PeacoQC"))
# }
# 
# fs <- flowCore::fsApply(flowSet, function(ff){
#   resQC_AI <- flow_auto_qc(fcsfiles = ff,
#                         folder_results = file.path(QC_dir, paste("flowAI", gsub(".fcs", " Folder", ff@description$GUID), sep = "/")),
#                         output = 1)
#   resQC <- PeacoQC(ff = ff,
#                    determine_good_cells = "all",
#                    channels = c(5:13),
#                    plot = TRUE,
#                    output_folder = file.path(QC_dir, "PeacoQC"))
#   ff <- ff[resQC$GoodCells, ]
#   return(ff)
# })

fs <- flowSet

fcs_colname <- colnames(fs)
antigen <- panel_md$antigen
marker_class <- panel_md$marker_class
panel <- as.data.frame(cbind(fcs_colname, antigen, marker_class))
fs[[1]]@description$`$CYT`
fs[[1]]@description$`$CYT` <- "FACS"

chs_of_interest <- colnames(fs)[5:13]
plot_aggregate(fs, channels = chs_of_interest, output_image = "FCSpreNorm.png")
normFlowSet <- warpSet(fs, stains = colnames(fs)[5])
plot_aggregate(normFlowSet, channels = chs_of_interest, output_image = "FCSpostNorm.png")

normFlowSet[[1]]@description$`$CYT`
panel$marker_class[11] <- "state"
sce <- CATALYST::prepData(normFlowSet, panel = panel, md = sample_md, transform = FALSE)
assay(sce, "exprs") <- assay(sce, "counts")
assays(sce)

CATALYST::pbMDS(sce, color_by = "condition", features = type_markers(sce), fun = "median", label_by = NULL)

p <- plotExprs(sce, features = NULL, color_by = "condition")
p$facet$params$ncol <- 4
p

n_events <- min(n_cells(sce))
n_events
n_cells(sce)
plotCounts(sce, group_by = "sample_id", color_by = "condition")

plotNRS(sce, features = type_markers(sce), color_by = "condition")

saveRDS(sce, file = "SCE_part1.rds")
