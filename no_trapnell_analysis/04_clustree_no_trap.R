
# Perform clustering at multiple resolutions and store in metadata
resolutions <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6,0.7, 0.8, 0.9, 1.0 , 1.1)

for (res in resolutions) {
  # Run FindNeighbors and FindClusters for each resolution
  seurat_integrated_genelist <- FindNeighbors(seurat_integrated_genelist, dims = 1:10, reduction = "pca")
  seurat_integrated_genelist <- FindClusters(seurat_integrated_genelist, resolution = res)
  
  # Ensure cluster identities are saved with a resolution-specific name
  seurat_integrated_genelist@meta.data[[paste0("RNA_snn_res.", res)]] <- Idents(seurat_integrated_genelist)
}

# View the metadata to confirm clusters for different resolutions
head(seurat_integrated_genelist@meta.data)

# Visualise the clustering using Clustree
library(clustree)
clustree(seurat_integrated_genelist@meta.data, prefix = "RNA_snn_res.")
