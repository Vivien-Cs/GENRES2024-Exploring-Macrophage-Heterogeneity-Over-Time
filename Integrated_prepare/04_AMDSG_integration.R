#Run AMDSG list integration#

seurat_integrated <- IntegrateLayers(
  object = seurat_merged, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

seurat_integrated <- FindNeighbors(seurat_integrated, features = gene_names)
seurat_integrated <- FindClusters(seurat_integrated, resolution = 0.6, cluster.name = "AMDSG_48hr_clusters")
seurat_integrated <- RunTSNE(seurat_integrated, features = gene_names)

DimPlot(seurat_integrated, label.box = T, label = T)
DimPlot(seurat_integrated, split.by = "condition_time")
#FeaturePlot(seurat_integrated, features = "Tnf", split.by = "condition_time")

