#Load packages
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

setwd("/mnt/scratch/users/vc680/single_cell/")

#save all files and folder within your working directory or working directory/your_custom_folder_name/
#update exp_path if you have a your_custom_folder_name within working directory 
# check if sub directory exists 

exp_path <- paste0(getwd(),"data/GSE158094_RAW/")

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
