# Set Idents to condition_time
Idents(seurat_integrated_genelist) <- "condition_time"

#Prepare the Seurat object for DE analysis
seurat_integrated_genelist <- PrepSCTFindMarkers(seurat_integrated_genelist)

# Find markers for all condition_time groups
markers_per_cluster <- FindAllMarkers(seurat_integrated_genelist, 
                                      only.pos = TRUE,  # Only positive markers
                                      min.pct = 0.25,   # Minimum expression in 25% of cells
                                      logfc.threshold = 0.25)  # Minimum log fold change

# View the top markers for each cluster
head(markers_per_cluster)

# Filter to include only significant DE genes (adjusted p-value < 0.05)
significant_markers <- subset(markers_per_cluster, p_val_adj < 0.05)

# View top significant markers
head(significant_markers)

# Specify the output directory for saving the markers
output_path <- "G:/My Drive/Viv_MACHS/R/markers_merged/"

# Save the markers to a CSV file
write.csv(significant_markers, file.path(output_path, "markers_per_condition_time.csv"), row.names = TRUE)



###################### loop condition vs condition ############################

# Get all unique condition_time groups
condition_time_groups <- unique(Idents(seurat_integrated_genelist))

# Loop through each pair of condition_time groups
for (i in seq_along(condition_time_groups)) {
  for (j in seq_along(condition_time_groups)) {
    if (i < j) {
      # Define the current pair of condition_time groups
      condition_time1 <- condition_time_groups[i]
      condition_time2 <- condition_time_groups[j]
      
      # Run DE analysis: condition_time1 vs. condition_time2
      de_genes <- FindMarkers(seurat_integrated_genelist, 
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


