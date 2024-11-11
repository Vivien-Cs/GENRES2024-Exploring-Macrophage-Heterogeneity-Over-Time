# **MACHS- Macrophage heterogeneity study**
## Project Description
This project investigates macrophage heterogeneity by analysing single-cell RNA sequencing
(scRNA-seq) data. It includes workflows for data integration, clustering, differential expression (DE) 
analysis, and visualisation. This project contains scripts for each analysis step.

## Project Structure

This project contains scripts for each analysis step.
Two versions of the workflow are available: 
`Integrated_prepare` (includes the Trapnell dataset) and `no_trap` (excludes the Trapnell dataset).

#Example of the no_trap workflow
Scripts are numbered for sequential execution:
- **`00_pkg_no_trap.R`**: Loads required packages.
- **`01_import_data_no_trap.R`**: Imports and prepares raw data files.
- **`02_preprocessing_no_trap.R`**: Preprocesses datasets, performing quality control and filtering.
- **`03_merge_and_integrate_no_trap.R`**: Merges and integrates Seurat objects.
- **`04_clustering_no_trap.R`** and **`04_clustree_no_trap.R`**: Clustering and visualisation of clustering stability.
- **`06_differential_expression_no_trap.R`**: Runs differential expression analysis.
- **`07_heatmap_no_trap.R`**: Generates heatmaps and performs DE analysis by cluster.
- **`08_line_plot_cluster.R`** and **`09_trend_analysis.R`**: Creates line plots and trend analysis for gene expression.

Output files are saved in the `output_integrated` directory, and figures are stored in `G:/My Drive/Viv_MACHS/R/`.


## Usage Instructions

### Running the Workflow
You can run the full workflow in two ways:
1. **Sequential Script Execution**: Run each script from `00` to `09` in numerical order.
2. **R Markdown Execution**: Run the entire analysis via the provided `.Rmd` file, which integrates
all steps for an automated workflow.

### Setting Up File Paths
Set `base_path` and `code_path` to match your directory structure
   
   
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
**`install.packages("readxl")`**,
**`install.packages("tidyr")`**,
**`install.packages("readxl")`**,
**`install.packages("circlize")`**,
**`install.packages("openxlsx")`**,
**`install.packages("clustree")`**


BioCManager is used to install and manage packages from the Bioconductor project for analysis of genomic data;  for the packages ComplexHeatmap,EnhancedVolcano and SingleCellExperiment you need BiocManager : 
**`if (!require("BiocManager", quietly = TRUE))`** 
**`install.packages("ComplexHeatmap")`**
**`BiocManager::install("EnhancedVolcano")`**
**`BiocManager::install("SingleCellExperiment")`**
**`BiocManager::install("slingshot")`**
**`BiocManager::install("MAST")`**


# Author and Contact

This project was created by **Vivien Csonka** under the supervision of **Shoumit Dey**. 
For questions or issues, please reach out via email: [shoumit.dey@york.ac.uk].
