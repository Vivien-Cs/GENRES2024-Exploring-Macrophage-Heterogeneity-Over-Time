
# Count the number of cells in each condition per cluster
# Extract metadata including clusters and conditions
cell_metadata <- seurat_integrated_genelist@meta.data

# Create a table to count the number of cells per cluster and condition
cell_counts <- cell_metadata %>%
  group_by(seurat_clusters, condition_time) %>%
  summarise(count = n()) %>%
  ungroup()

# Calculate proportions of cells per condition in each cluster
# Total number of cells in each condition
total_cells_per_condition <- cell_metadata %>%
  group_by(condition_time) %>%
  summarise(total_cells = n()) 

# Merge with the cell counts table to calculate proportions
cell_proportions <- cell_counts %>%
  left_join(total_cells_per_condition, by = "condition_time") %>%
  mutate(proportion = count / total_cells)

#  Reshape the data into a matrix format for heatmap plotting
# Convert the data to a wide format, with clusters as rows and conditions as columns
heatmap_data <- cell_proportions %>%
  select(seurat_clusters, condition_time, proportion) %>%
  tidyr::pivot_wider(names_from = condition_time, values_from = proportion, values_fill = 0)

# Reorder the condition_time manually to match the desired order
time_order <- c("Media_48_Li_2019_zero","Media_M0_Hagai_2018_zero",
                "LPS_M1_Hagai_2018_two", "LPS_M1_Hagai_2018_four", "LPS_M1_Hagai_2018_six", 
                "M1_LPS_IFNG_48_Li_2019_fourty8")

heatmap_data <- heatmap_data %>%
  select(seurat_clusters, all_of(time_order)) # Reordering the columns based on the sample order

# Convert to matrix for heatmap plotting
heatmap_matrix <- as.matrix(heatmap_data[, -1])  # Remove cluster column 0s

# create row labels based on the number of clusters
num_clusters <- nrow(heatmap_matrix)
row_labels <- paste0("Cluster ", 1:num_clusters)  # Automatically generate cluster labels based on row count

# Plot the heatmap with a border and properly ordered cluster labels
col_fun <- colorRamp2(c(0, max(heatmap_matrix, na.rm = TRUE)), c("white", "blue"))

ploth <- ComplexHeatmap::Heatmap(heatmap_matrix,
                        name = "Proportion of Cells",
                        col = col_fun,
                        cluster_rows = TRUE,   
                        cluster_columns = FALSE, 
                        show_row_names = TRUE, 
                        show_column_names = TRUE, 
                        row_title = "Clusters",  
                        column_title = "Samples", 
                        row_names_side = "left",
                        row_names_gp = gpar(fontsize = 10),
                        column_names_gp = gpar(fontsize = 10), 
                        row_labels = row_labels, 
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.rect(x = x, y = y, width = width, height = height,
                                    gp = gpar(col = "black", fill = NA)) 
                        }
)

# Use png() to specify the file and dimensions
png("G:/My Drive/Viv_MACHS/R/heatmap_no_Trap_res_0.7.png", width = 10, height = 8, units = "in", res = 300)

# Print the heatmap to save it in the specified device
draw(ploth)

# Close the device
dev.off()

#################### Find DE per cluster #######################
#  Set Idents to seurat_clusters to find DE genes for each cluster
Idents(seurat_integrated_genelist) <- "seurat_clusters"
# Rename clusters from 0-based to 1-based numbering
current_clusters <- levels(seurat_integrated_genelist)
new_cluster_labels <- as.character(as.numeric(current_clusters) + 1)

#  Set new cluster identities
seurat_integrated_genelist <- RenameIdents(seurat_integrated_genelist, 
                                           setNames(new_cluster_labels, current_clusters))

# Verify the cluster renaming
levels(seurat_integrated_genelist)  


seurat_integrated_genelist <- PrepSCTFindMarkers(seurat_integrated_genelist)

# run the DE analysis 
markers_per_cluster <- FindAllMarkers(seurat_integrated_genelist, 
                                      only.pos = TRUE,  
                                      min.pct = 0.25,   
                                      logfc.threshold = 0.25)


#  Filter to include only significant DE genes (adjusted p-value < 0.05)
significant_markers <- markers_per_cluster %>%
  filter(p_val_adj < 0.05)

# View the top markers for each cluster
head(significant_markers)

#  Split the significant markers by cluster
markers_by_cluster <- split(significant_markers, significant_markers$cluster)

#Save the markers for each cluster to separate sheets in an Excel file
output_path <- "G:/My Drive/Viv_MACHS/R/markers_per_cluster_no_trapnell_0.7.xlsx"
wb <- createWorkbook()

for (cluster_name in names(markers_by_cluster)) {
  # Add a new sheet for each cluster
  addWorksheet(wb, paste0("Cluster_", cluster_name))
  
  # Write the data for each cluster to the corresponding sheet
  writeData(wb, sheet = paste0("Cluster_", cluster_name), markers_by_cluster[[cluster_name]])
}

# Save the workbook
saveWorkbook(wb, output_path, overwrite = TRUE)


FeaturePlot(seurat_integrated_genelist, features = "Tnf", 
            reduction = "umap", 
            label = TRUE,
            split.by = "condition_time")


