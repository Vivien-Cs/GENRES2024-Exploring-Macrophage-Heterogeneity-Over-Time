## slingshot
seurat_merged <- readRDS(paste0("preprocessed_merged_dims4_",res,".Rds"))
# interpet trajectories
#Load packages

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