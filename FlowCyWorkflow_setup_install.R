rm(list = ls())

install_packages <- FALSE

if(!require(devtools)){
  BiocManager::install("devtools")
}

if(install_packages){
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)
  devtools::install_github("VPetukhov/ggrastr")
  devtools::install_github("immunogenomics/harmony")
  BiocManager::install("cydar")
  BiocManager::install("scuttle")
  BiocManager::install("GenomeInfoDb")
  devtools::install_github("casanova-lab/iMUBAC")
  devtools::install_github("masato-ogishi/plotUtility")
  BiocManager::install("flowAI")
  remotes::install_github("saeyslab/CytoNorm")
  BiocManager::install("PeacoQC")
  BiocManager::install("ggpubr")
}

if(install_packages){
  if(!require('rstudioapi')) {
    install.packages('rstudioapi')
  }
  if(!require('devtools')){
    install.packages("devtools")
  }
  if(!require('flowCore')){
    BiocManager::install("flowCore")
  }
  if(!require('cytofCore')){
    devtools::install_github("nolanlab/cytofCore")
  }
  # if(!require('JinmiaoChenLab/cytofkit')){
  #   #do not update hexbin
  #   remotes::install_github("JinmiaoChenLab/cytofkit")
  # }
  # if(!require("CytoML")){
  #   BiocManager::install("CytoML")
  # }
  if(!require('FlowSOM')){
    BiocManager::install("FlowSOM")
  }
  if(!require('cluster')){
    install.packages("cluster")
  }
  if(!require('Rtsne')){
    install.packages("Rtsne")
  }
  if(!require('ggplot2')){
    install.packages("ggplot2")
  }
  if(!require('dplyr')){
    install.packages("dplyr")
  }
  if(!require('ggthemes')){
    install.packages("ggthemes")
  }
  if(!require('RColorBrewer')){
    install.packages('RColorBrewer')
  }
  if(!require("uwot")){
    install.packages("uwot")
  }
  if(!require("CATALYST"))
    BiocManager::install("CATALYST")
  if(!require("diffcyt"))
    BiocManager::install("diffcyt")
  if(!require("stringr"))
    BiocManager::install("stringr")
  
  if(!require("JinmiaoChenLab/Rphenograph")){
    remotes::install_github("JinmiaoChenLab/Rphenograph")
  }
  BiocManager::install("Rphenograph")
  if(!require("scran"))
    BiocManager::install("scran")
  if(!require("scater"))
    BiocManager::install("scater")
  if(!require("ggcyto"))
    BiocManager::install("ggcyto")
  if(!require("SingleCellExperiment"))
    BiocManager::install("SingleCellExperiment")
  if(!require("Rphenograph"))
    BiocManager::install("Rphenograph")
  if(!require("flowWorkspace"))
    BiocManager::install("flowWorkspace")
  if(!require("flowVS"))
    install.packages(file.choose(), repos = NULL, type = "source")
  if(!require("flowStats"))
    BiocManager::install("flowStats")
}

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
options(java.parameters="-Xmx60G")
library(tidyverse)
library(data.table)
library(scuttle)
library(iMUBAC)
