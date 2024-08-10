
# Perform DE analysis for each cluster ##
cluster_markers <- list()
for (i in 0:5) {  # Assuming clusters are labeled from 0 to 5
  cluster_markers[[paste0("Cluster_", i)]] <- FindMarkers(seurat_merged, ident.1 = i, min.pct = 0.25, logfc.threshold = 0.25)
}

# Print results for each cluster
for (i in 0:5) {
  cat(paste0("Markers for Cluster ", i, ":\n"))
  print(cluster_markers[[paste0("Cluster_", i)]])
  cat("\n")
}

library(openxlsx)
# Write results for each cluster to an Excel file
output_directory <- "/mnt/scratch/users/vc680/single_cell/R/markers_merged/"
for (i in 0:5) {
  markers_cluster <- paste0(output_directory, "/Cluster_", i, "_markers.xlsx")
  write.xlsx(cluster_markers[[paste0("Cluster_", i)]], file = markers_cluster, row.names = TRUE)
  cat(paste0("Markers for Cluster ", i, " written to ", markers_cluster, "\n"))
}








###################################################################################################


#  Create a new identity that combines condition and time point
seurat_merged$condition_time <- paste(seurat_merged$condition, seurat_merged$time_pt, sep = "_")
Idents(seurat_merged) <- "condition_time"

#  Define the groups for comparison
M0_group <- c("M0_zero", "M0_six", "M0_twenty4")
M1_group <- c("M1_six", "M1_twenty4","M1_fourty8")

#  Find markers between M0 and M1 groups
de_genes <- FindMarkers(seurat_merged, 
                        ident.1 = M0_group, 
                        ident.2 = M1_group,
                        
                        only.pos = TRUE,   # Include both up and down-regulated genes if only.pos=false
                        min.pct = 0.25, 
                        logfc.threshold = 0.25)

de_genes <- subset(de_genes, subset = p_val_adj < 0.05)


# Save the markers
write.csv(de_genes, file.path(output_integrated, "markers_merged", paste0("DE_genes_M0_vs_M1_allTimepoints_dims", dims, "_res", res, ".csv")), row.names = TRUE)

# Read the markers from the CSV file
all_markers <- read.csv(file.path(output_integrated, "markers_merged", paste0("DE_genes_M0_vs_M1_allTimepoints_dims", dims, "_res", res, ".csv")), row.names = 1)

# Select top markers
top20 <- all_markers %>% 
  arrange(desc(abs(avg_log2FC))) %>% 
  slice_head(n = 20)

# Create a heatmap 
DoHeatmap(seurat_merged, features = rownames(top20), group.by = "condition_time") + 
  NoLegend()

# # Save the heatmap
ggsave(file.path(output_integrated, "markers_merged", "DE_genes_M0_vs_M1_allTimepoints_heatmap.pdf"), width = 15, height = 20)
