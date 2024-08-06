# **Stem Cell RNA-Seq Analysis**

# Description


# Packages
###Packages required to run data:
**`install.packages("Seurat")`**,
**`install.packages("SeuratObject")`**,
**`install.packages("cowplot")`** 
**`install.packages("dplyr")`**,
**`install.packages("ggplot2")`**, 
**`install.packages("stringr")`**, 
**`install.packages("ggpubr")`**,
**`install.packages ("reshape2")`**,
**`install.packages("colorRamp2")`**, 
**`install.packages("kableExtra")`**,
**`install.packages("bookdown")`**,
**`install.packages("ggridges")`**,
**`install.packages("sqldf")`**,
**`install.packages("readxl")`**

BioCManager is used to install and manage packages from the Bioconductor project for analysis of genomic data;  for the packages ComplexHeatmap,EnhancedVolcano and SingleCellExperiment you need BiocManager : 
**`if (!require("BiocManager", quietly = TRUE))`** 
**`install.packages("ComplexHeatmap")`**
**`BiocManager::install("EnhancedVolcano")`**
**`BiocManager::install("SingleCellExperiment")`**
**`BiocManager::install("slingshot")`**


