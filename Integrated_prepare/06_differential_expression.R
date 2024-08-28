# 
# # Perform DE analysis for each cluster ##
# cluster_markers <- list()
# for (i in 0:5) {  # Assuming clusters are labeled from 0 to 5
#   cluster_markers[[paste0("Cluster_", i)]] <- FindMarkers(seurat_merged, ident.1 = i, min.pct = 0.25, logfc.threshold = 0.25)
# }
# 
# # Print results for each cluster
# for (i in 0:5) {
#   cat(paste0("Markers for Cluster ", i, ":\n"))
#   print(cluster_markers[[paste0("Cluster_", i)]])
#   cat("\n")
# }
# 
# library(openxlsx)
# # Write results for each cluster to an Excel file
# output_directory <- "~/drive-download-20240617T134729Z-001/R/markers_merged/"
# for (i in 0:5) {
#   markers_cluster <- paste0(output_directory, "/Cluster_", i, "_markers.xlsx")
#   write.xlsx(cluster_markers[[paste0("Cluster_", i)]], file = markers_cluster, row.names = TRUE)
#   cat(paste0("Markers for Cluster ", i, " written to ", markers_cluster, "\n"))
# }
# 

###################################################################################################
#RUN DE MAST#
# Assume seurat_integrated is your integrated Seurat object
Idents(seurat_integrated) <- "condition_time"

#Prepare the Seurat object for DE analysis
seurat_integrated <- PrepSCTFindMarkers(seurat_integrated)

# Define the groups for comparison
M0_group <- c("M0_zero", "M0_six", "M0_twenty4")
M1_group <- c("M1_six", "M1_twenty4", "M1_fourty8")

# Find markers between M0 and M1 groups using MAST on the integrated dataset
de_genes <- FindMarkers(seurat_integrated, 
                        ident.1 = M0_group, 
                        ident.2 = M1_group,
                        test.use = "MAST",    
                        only.pos = TRUE,      
                        min.pct = 0.25, 
                        logfc.threshold = 0.25)

# Filter, save, and visualize as before
de_genes <- subset(de_genes, subset = p_val_adj < 0.05)

# Specify the output directory
output_path <- "G:/My Drive/Genres2024/output/"

# Save the DE genes as a CSV file named "markers_integrated_amdsg.csv"
write.csv(de_genes, file.path(output_path, "markers_integrated_condition_timeamdsg.csv"), row.names = TRUE)

all_markers <- read.csv(file.path(output_integrated, "markers_merged", paste0("DE_genes_M0_vs_M1_allTimepoints_dims", dims, "_res", res, ".csv")), row.names = 1)
top20 <- all_markers %>% arrange(desc(abs(avg_log2FC))) %>% slice_head(n = 20)
DoHeatmap(seurat_integrated, features = rownames(top20), group.by = "condition_time") + NoLegend()
ggsave(file.path(output_integrated, "markers_merged", "DE_genes_M0_vs_M1_allTimepoints_heatmap.pdf"), width = 15, height = 20)



############ Loop through condition_time vs condition_time ###########
# Set the identity class to condition_time
Idents(seurat_integrated) <- "condition_time"

# Prepare the Seurat object for DE analysis
seurat_integrated <- PrepSCTFindMarkers(seurat_integrated)

# Get all unique condition_time groups
condition_time_groups <- unique(Idents(seurat_integrated))

# Specify the output directory
output_path <- "G:/My Drive/Genres2024/output/DE_analysis/"

# Loop through each pair of condition_time groups
for (i in seq_along(condition_time_groups)) {
  for (j in seq_along(condition_time_groups)) {
    if (i < j) {
      # Define the current pair of condition_time groups
      condition_time1 <- condition_time_groups[i]
      condition_time2 <- condition_time_groups[j]
      
      # Run DE analysis: condition_time1 vs. condition_time2
      de_genes <- FindMarkers(seurat_integrated, 
                              ident.1 = condition_time1, 
                              ident.2 = condition_time2, 
                              test.use = "MAST",    
                              only.pos = TRUE,      
                              min.pct = 0.25, 
                              logfc.threshold = 0.25)
      
      # Filter to include only significant genes
      de_genes <- subset(de_genes, subset = p_val_adj < 0.05)
      
      # Save the DE genes as a CSV file
      output_filename <- paste0("DE_genes_", condition_time1, "_vs_", condition_time2, ".csv")
      write.csv(de_genes, file.path(output_path, output_filename), row.names = TRUE)
    }
  }
}


########## now loop through condition vs condition ###############

# Set the identity class to condition
Idents(seurat_integrated) <- "condition"

# Get all unique time points
time_points <- unique(seurat_integrated$time_pt)

# Specify the output directory
output_path <- "G:/My Drive/Genres2024/output/DE_analysis/Condition_vs_Condition/"

# Loop through each time point and compare conditions within that time point
for (time in time_points) {
  
  # Subset the Seurat object for the current time point
  seurat_time_subset <- subset(seurat_integrated, subset = time_pt == time)
  
  # Prepare the subsetted Seurat object for DE analysis
  seurat_time_subset <- PrepSCTFindMarkers(seurat_time_subset)
  
  # Get the unique conditions within the current time point
  conditions <- unique(Idents(seurat_time_subset))
  
  # Perform pairwise comparisons between conditions at this time point
  for (i in seq_along(conditions)) {
    for (j in seq_along(conditions)) {
      if (i < j) {
        # Define the current pair of conditions
        condition1 <- conditions[i]
        condition2 <- conditions[j]
        
        # Run DE analysis: condition1 vs. condition2 at the current time point
        de_genes <- FindMarkers(seurat_time_subset, 
                                ident.1 = condition1, 
                                ident.2 = condition2, 
                                test.use = "MAST",    
                                only.pos = TRUE,      
                                min.pct = 0.25, 
                                logfc.threshold = 0.25)
        
        # Filter to include only significant genes
        de_genes <- subset(de_genes, subset = p_val_adj < 0.05)
        
        # Save the DE genes as a CSV file
        output_filename <- paste0("DE_genes_", condition1, "_vs_", condition2, "_at_", time, ".csv")
        write.csv(de_genes, file.path(output_path, output_filename), row.names = TRUE)
      }
    }
  }
}


############### now loop through time point vs time point #########################
# Set the identity class to time_pt
Idents(seurat_integrated) <- "time_pt"

# Prepare the Seurat object for DE analysis
seurat_integrated <- PrepSCTFindMarkers(seurat_integrated)

# Get all unique conditions
conditions <- unique(seurat_integrated$condition)

# Specify the output directory
output_path <- "G:/My Drive/Genres2024/output/DE_analysis/Time_vs_Time/"

# Loop through each condition and compare time points within that condition
for (condition in conditions) {
  
  # Subset the Seurat object for the current condition
  seurat_condition_subset <- subset(seurat_integrated, subset = condition == condition)
  
  # Get the unique time points within the current condition
  time_points <- unique(Idents(seurat_condition_subset))
  
  # Perform pairwise comparisons between time points in this condition
  for (i in seq_along(time_points)) {
    for (j in seq_along(time_points)) {
      if (i < j) {
        # Define the current pair of time points
        time_point1 <- time_points[i]
        time_point2 <- time_points[j]
        
        # Run DE analysis: time_point1 vs. time_point2 within the current condition
        de_genes <- FindMarkers(seurat_condition_subset, 
                                ident.1 = time_point1, 
                                ident.2 = time_point2, 
                                test.use = "MAST",    
                                only.pos = TRUE,      
                                min.pct = 0.25, 
                                logfc.threshold = 0.25)
        
        # Filter to include only significant genes
        de_genes <- subset(de_genes, subset = p_val_adj < 0.05)
        
        # Save the DE genes as a CSV file
        output_filename <- paste0("DE_genes_", time_point1, "_vs_", time_point2, "_in_", condition, ".csv")
        write.csv(de_genes, file.path(output_path, output_filename), row.names = TRUE)
      }
    }
  }
}
