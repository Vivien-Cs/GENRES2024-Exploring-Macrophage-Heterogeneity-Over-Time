### Packages
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
base_path <- "G:/My Drive/Genres2024"

# Define the output directory
output_integrated <- "G:/My Drive/Genres2024/output/"
if (!dir.exists(output_integrated)) {
  dir.create(output_integrated, recursive = TRUE)
}

# Function to process and analyse a single dataset without integration
process_single_dataset <- function(lstSamples_2_4_6h, lstConditions_2_4_6h, lstTime_2_4_6hr, data_dir_base, data_set_name = "Hagai_2018_amdsg_psg") {
  #create a list to hold all the objects within a dataset
  lstSeuratObj <- list()
      for (sample in 1:length(lstSamples_2_4_6h)) {
          # Construct the correct data directory path for the current sample
          data_dir <- file.path(data_dir_base, lstConditions_2_4_6h[[sample]], lstSamples_2_4_6h[[sample]])
          print(data_dir)
          
          # Load the data for the current sample
          seuratObj <- Read10X(data.dir = data_dir)
          seuratObj <- CreateSeuratObject(counts = seuratObj, project = lstSamples_2_4_6h[[sample]], min.cells = 3, min.features = 200)
          seuratObj$condition <- lstConditions_2_4_6h[[sample]]
          seuratObj$time_pt <- lstTime_2_4_6hr[[sample]]
          
          # Preprocessing steps
          seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = "^mt-")
          mito_genes <- grep("^mt-", rownames(seuratObj), value = TRUE)
          print(paste("Mitochondrial genes detected:", mito_genes))
          
          # Visualization and quality control
          p1 <- VlnPlot(seuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
          p2 <- FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
          print(p1 + p2)
          
          #subset
          seuratObj <- subset(seuratObj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
          
          #merge by condition and time
          seuratObj$condition_time <- paste0(seuratObj$condition,"_",seuratObj$time_pt)
          #add to list
          lstSeuratObj[[lstSamples_2_4_6h[[sample]]]] <- seuratObj
          print(paste0("Saved list item:",lstSamples_2_4_6h[[sample]]))
      }
  #merging all samples within the dataset
  seurat_merged <- merge(
    lstSeuratObj[[lstSamples_2_4_6h[[1]]]], 
    y = lstSeuratObj[2:length(lstSeuratObj)]
  )
  
  dims=30
  #############################################
  #################################################
  # Normalisation with SCTransform
  #################################################
  ##################################################
 #seurat_merged <- SCTransform(seurat_merged, vst.flavor = "v2")
  #############################################
  #################################################
  # Normalisation with older workflow
  #################################################
  seurat_merged<-NormalizeData(seurat_merged)
  seurat_merged<-ScaleData(seurat_merged)
  
  #PSG list
  PSG <- read_excel("G:/My Drive/Genres2024/R/PSG_48hr.xlsx")
  psg_gene <- PSG$Gene
  # #AMDSG
  AMDSG <- read_excel("G:/My Drive/Genres2024/R/AMDSG_48hr.xlsx")
  amdsg_gene <- AMDSG$Gene
  # 
  #combined
  gene_names <- unique(c(psg_gene,amdsg_gene))
  
  # Dimensionality reduction with PCA and UMAP
  seurat_merged <- RunPCA(seurat_merged, features = gene_names)
  seurat_merged <- RunUMAP(seurat_merged, dims = 1:dims)
      
  # Clustering
  seurat_merged <- FindNeighbors(seurat_merged, dims = 1:dims)
  seurat_merged <- FindClusters(seurat_merged, resolution = 0.5)
  
  # Generate and save UMAP plots
  all_cell <- DimPlot(seurat_merged, label.box = TRUE, label = TRUE)
  ggsave(paste0(output_integrated, "all_cells_", data_set_name, ".pdf"), plot = all_cell, width = 10, height = 8)
  
  split_plot <- DimPlot(seurat_merged, split.by = "condition_time", ncol = 5)
  ggsave(paste0(output_integrated, "split_by_condition_time_", data_set_name, ".pdf"), plot = split_plot, width = 15, height = 20)
  
  split_plot <- DimPlot(seurat_merged, split.by = "condition", ncol = 5)
  ggsave(paste0(output_integrated, "split_by_condition_", data_set_name, ".pdf"), plot = split_plot, width = 15, height = 20)
  
  split_plot <- DimPlot(seurat_merged, split.by = "orig.ident", ncol = 5)
  ggsave(paste0(output_integrated, "split_by_origin_ident_", data_set_name, ".pdf"), plot = split_plot, width = 15, height = 20)
}

# Example: Process a single dataset (48hr)
lstSamples_2_4_6h <- c("LPS_M1_2h_rep1", "LPS_M1_4h_rep1", "LPS_M1_6h_rep1",
                       "LPS_M1_2h_rep2", "LPS_M1_4h_rep2", "LPS_M1_6h_rep2",
                       "LPS_M1_2h_rep3", "LPS_M1_4h_rep3", "LPS_M1_6h_rep3",
                       "Med_M0_0h_rep1", "Med_M0_0h_rep2", "Med_M0_0h_rep3")

lstConditions_2_4_6h <- c("LPS_M1_Hagai_2018", "LPS_M1_Hagai_2018", "LPS_M1_Hagai_2018",
                          "LPS_M1_Hagai_2018", "LPS_M1_Hagai_2018", "LPS_M1_Hagai_2018",
                          "LPS_M1_Hagai_2018", "LPS_M1_Hagai_2018", "LPS_M1_Hagai_2018",
                          "Media_M0_Hagai_2018", "Media_M0_Hagai_2018", "Media_M0_Hagai_2018")
lstTime_2_4_6hr <- c("two", "four", "six","two", "four", "six","two", "four", "six","zero", "zero", "zero")

data_dir_base <- paste0(base_path,"/data_raw/E_MTAB_6754_2_4_6hr_RAW")
process_single_dataset(lstSamples_2_4_6h, lstConditions_2_4_6h,lstTime_2_4_6hr, data_dir_base, data_set_name = "Hagai_2018_amdsg_psg")
