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
library(ComplexHeatmap)
library(colorRamp2)
library(ggridges)
library(slingshot)
library(SingleCellExperiment)

# Set the working directory to your project folder.
setwd("G:/My Drive/Genres2024")

# Define the base path for the project.
base_path<- "G:/My Drive/Genres2024"
#save all files and folder within your working directory or working directory/your_custom_folder_name/
#update exp_path if you have a your_custom_folder_name within working directory 
# check if sub directory exists 

## Ensure that all files and folders are saved within your working directory or a custom sub-folder.
## If using a custom sub-folder, update `exp_path` accordingly.
exp_path <- paste0(getwd(),"/data_raw/GSE158094_RAW/")


if (!file.exists(paste0(getwd(), "/R"))){
  dir.create(paste0(getwd(), "/R"))
}

output_integrated <- paste0(getwd(), "/R/")

#set pca dimensions to use and 
#resolution for cluster identification
dims=10
res=0.3

#create sub-directories
dir.create(paste0(output_integrated, "group_comparisons"))
dir.create(paste0(output_integrated, "markers"))
dir.create(paste0(output_integrated, "markers_merged"))
dir.create(paste0(output_integrated, "phenotyping"))

