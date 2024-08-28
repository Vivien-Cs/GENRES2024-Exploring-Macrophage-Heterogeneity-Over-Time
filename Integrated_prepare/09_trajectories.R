
# Ensure all M0_0h samples are labeled consistently
seurat_merged$condition_time[seurat_merged$condition_time %in% c("M0_0h", "M0_zero", "M0_0h_48")] <- "M0_0h"

# Extract PCA and cluster labels directly from Seurat object
rd <- Embeddings(seurat_merged, "pca")
cl <- seurat_merged$seurat_clusters

# Ensure correct format
rd <- as.matrix(rd)
cl <- as.factor(cl)

# Identify all cells from M0_0h condition
m0_0hr_cells <- WhichCells(seurat_merged, expression = condition_time == "M0_0h")

# Get all clusters containing M0_0h cells
m0_0hr_clusters <- unique(seurat_merged$seurat_clusters[m0_0hr_cells])

# Convert to character vector
start_clusters <- as.character(m0_0hr_clusters)

# Print the starting clusters for verification
print(paste("Starting clusters:", paste(start_clusters, collapse = ", ")))

# Run Slingshot
set.seed(123)  
slingshot_obj <- slingshot(rd, clusterLabels = cl, start.clus = start_clusters)

# Add pseudotime to Seurat object
seurat_merged$pseudotime <- slingPseudotime(slingshot_obj)

# Run t-SNE
if (!"tsne" %in% Reductions(seurat_merged)) {
  seurat_merged <- RunTSNE(seurat_merged, dims = 1:18)
}

# Extract t-SNE coordinates
tsne_coords <- Embeddings(seurat_merged, "tsne")

# Create a data frame for ggplot
tsne_plot_df <- data.frame(tsne_coords, 
                           cluster = seurat_merged$seurat_clusters, 
                           pseudotime = seurat_merged$pseudotime,
                           condition_time = seurat_merged$condition_time)

# Plotting with ggplot2
p <- ggplot(tsne_plot_df, aes(x = tSNE_1, y = tSNE_2, color = pseudotime)) +
  geom_point(size = 1.5) +
  scale_color_gradient(low = "yellow", high = "blue") +
  theme_minimal() +
  labs(title = "t-SNE with Slingshot Pseudotime", x = "t-SNE1", y = "t-SNE2") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  theme(legend.position = "right")

# Add text labels for M0_0h cells
m0_0h_cells_df <- tsne_plot_df[tsne_plot_df$condition_time == "M0_0h", ]
p <- p + geom_text(data = m0_0h_cells_df, aes(label = "M0_0h"), color = "red", size = 3, check_overlap = TRUE)

# Save the plot
ggsave("G:/My Drive/Genres2024/output/tsne_pseudotime_plot.png", plot = p, width = 10, height = 8, dpi = 300)

# Display the plot
print(p)




# 
# # Convert Seurat to SingleCellExperiment
#   sce <- as.SingleCellExperiment(seurat_merged)
# 
# # Extract PCA and cluster labels
# rd <- reducedDim(sce, "PCA")   
# cl <- sce$seurat_clusters
# 
# # Ensure correct format (fixed)
# rd <- as.matrix(rd)
# cl <- as.factor(cl)
# 
# # Identify the starting cluster for M0_0hr
# m0_0hr_cells <- WhichCells(seurat_merged, expression = condition_time == "M0_0h")
# m0_0hr_cluster <- unique(seurat_merged$seurat_clusters[m0_0hr_cells])
# start_cluster <- as.character(m0_0hr_cluster[1])  # Using the first matching cluster
# 
# # Run Slingshot
# set.seed(123)  
# slingshot_obj <- slingshot(rd, clusterLabels = cl, start.clus = start_cluster)
# 
# # Add pseudotime to Seurat object
# seurat_merged$pseudotime <- slingPseudotime(slingshot_obj)
# 
# # Run t-SNE
# seurat_merged <- RunTSNE(seurat_merged, dims = 1:18)
# 
# # Extract t-SNE coordinates
# tsne_coords <- Embeddings(seurat_merged, "tsne")
# 
# # Create a data frame for ggplot
# tsne_plot_df <- data.frame(tsne_coords, cluster = seurat_merged$seurat_clusters, pseudotime = seurat_merged$pseudotime)
# 
# # Extract curves for plotting
# curve_df <- data.frame()
# for (i in seq_along(slingCurves(slingshot_obj))) {
#   curve_i <- as.data.frame(slingCurves(slingshot_obj)[[i]]$s)
#   curve_i$curve_id <- i
#   curve_df <- rbind(curve_df, curve_i)
# }
# 
# # Plotting with ggplot2
# p <- ggplot(tsne_plot_df, aes(x = tSNE_1, y = tSNE_2, color = pseudotime)) +
#   geom_point(size = 1.5) +
#   scale_color_gradient(low = "yellow", high = "blue") +
#   geom_path(data = curve_df, aes(x = Dim1, y = Dim2, group = curve_id), color = "black", linewidth = 1) +
#   theme_minimal() +
#   labs(title = "t-SNE with Slingshot Pseudotime", x = "t-SNE1", y = "t-SNE2") +
#   theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
#   theme(legend.position = "right")
# 
# # Display the plot
# print(p)




# # Add pseudotime to Seurat object
# seurat_merged$pseudotime <- slingPseudotime(slingshot_obj)
# # Plot trajectory
# plot(rd[,1:2], col = cl, pch = 16, asp = 1)
# 
# # Extract and plot the lineages
# lineages <- slingLineages(slingshot_obj)
# curves <- slingCurves(slingshot_obj)
# 
# for (i in seq_along(curves)) {
#   lines(curves[[i]]$s[curves[[i]]$ord, 1:2], col = "black", lwd = 2)
# }
# 
# # Add legend
# legend("topright", legend = levels(cl), col = 1:length(levels(cl)), pch = 16)
# 
# 
# #pseudotime heatmap
# # Assign pseudotime based on the trajectory
# seurat_merged$pseudotime <- slingPseudotime(slingshot_obj)
# 
# # Use pseudotime in further analyses
# expression_data <- GetAssayData(seurat_merged, assay = "RNA", slot = "data")
# 
# # Check if the data is not matrix
# if (!is.matrix(expression_data)) {
#   # Convert to matrix 
#   expression_data <- as.matrix(expression_data)
# }
# 
# 
# # Order the cells by pseudotime and visualise gene expression changes
# ordered_cells <- order(seurat_merged$pseudotime, na.last = NA)
# 
# # Subset the expression data matrix by ordered cells
# expression_data_ordered <- expression_data[, ordered_cells]
# 
# 
# # Define color scale
# col_fun <- colorRamp2(c(min(expression_data_ordered), median(expression_data_ordered), max(expression_data_ordered)), c("blue", "white", "red"))
# 
# # Create and visualise a heatmap based on pseudotime ordering
# heatmap <- Heatmap(expression_data_ordered, name = "Expression", col = col_fun, use_raster = TRUE)
# 
# # Save the heatmap
# ggsave("heatmap_results/total_heatmap.png", plot = grid.grabExpr(draw(heatmap, heatmap_legend_side = "bottom")), device = "png", width = 10, height = 8, units = "in", dpi = 300)
# 
# 
# draw(heatmap)