
seurat_merged <- readRDS(paste0("preprocessed_merged_dims4_",res,".Rds"))
# interpet trajectories

# Convert Seurat to SingleCellExperiment
sce <- as.SingleCellExperiment(seurat_merged)

# Extract PCA and cluster labels
rd <- reducedDim(sce, "PCA")   
cl <- sce$seurat_clusters

# Ensure correct format (fixed)
rd <- as.matrix(rd)
cl <- as.factor(cl)

# Run Slingshot
set.seed(123)  
slingshot_obj <- slingshot(rd, clusterLabels = cl, start.clus = "0")

# Add pseudotime to Seurat object
seurat_merged$pseudotime <- slingPseudotime(slingshot_obj)
# Plot trajectory
plot(rd[,1:2], col = cl, pch = 16, asp = 1)

# Extract and plot the lineages
lineages <- slingLineages(slingshot_obj)
curves <- slingCurves(slingshot_obj)

for (i in seq_along(curves)) {
  lines(curves[[i]]$s[curves[[i]]$ord, 1:2], col = "black", lwd = 2)
}

# Add legend
legend("topright", legend = levels(cl), col = 1:length(levels(cl)), pch = 16)


#pseudotime heatmap
# Assign pseudotime based on the trajectory
seurat_merged$pseudotime <- slingPseudotime(slingshot_obj)

# Use pseudotime in further analyses
expression_data <- GetAssayData(seurat_merged, assay = "RNA", slot = "data")

# Check if the data is not matrix
if (!is.matrix(expression_data)) {
  # Convert to matrix 
  expression_data <- as.matrix(expression_data)
}


# Order the cells by pseudotime and visualise gene expression changes
ordered_cells <- order(seurat_merged$pseudotime, na.last = NA)

# Subset the expression data matrix by ordered cells
expression_data_ordered <- expression_data[, ordered_cells]


# Define color scale
col_fun <- colorRamp2(c(min(expression_data_ordered), median(expression_data_ordered), max(expression_data_ordered)), c("blue", "white", "red"))

# Create and visualise a heatmap based on pseudotime ordering
heatmap <- Heatmap(expression_data_ordered, name = "Expression", col = col_fun, use_raster = TRUE)

# Save the heatmap
ggsave("heatmap_results/total_heatmap.png", plot = grid.grabExpr(draw(heatmap, heatmap_legend_side = "bottom")), device = "png", width = 10, height = 8, units = "in", dpi = 300)


draw(heatmap)