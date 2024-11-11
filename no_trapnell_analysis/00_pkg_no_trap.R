## Check the README.md and sessioninfo.md files for additional information and instructions.
## Ensure that all required files are placed in the correct directories as outlined below.

# Load necessary packages for the analysis.
library(Seurat)
library(SeuratObject)
library(cowplot)
library(dplyr)
library(ggplot2)
library(stringr)
library(sqldf)
library(ggpubr)
library(reshape2)
library(EnhancedVolcano)
library(readxl)
library(slingshot)
library(viridis)
library(monocle3)
library(slingshot)
library(ComplexHeatmap)
library(circlize)
library(tidyr)
library(openxlsx)
library(harmony)
# library(colorRamp2)
# library(ggridges)
# library(SingleCellExperiment)

if (!file.exists(paste0(base_path, "/R"))){
  dir.create(paste0(base_path, "/R"))
}

output_integrated <- paste0(base_path, "R/")

#set pca dimensions to use and 
#resolution for cluster identification
dims=10
res=0.3

#create sub-directories
dir.create(paste0(output_integrated, "group_comparisons"))
dir.create(paste0(output_integrated, "markers"))
dir.create(paste0(output_integrated, "markers_merged"))
dir.create(paste0(output_integrated, "phenotyping"))
