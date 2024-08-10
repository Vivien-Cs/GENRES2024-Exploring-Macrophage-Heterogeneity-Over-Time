##integrate using entire dataset##

seurat_integrated <- IntegrateLayers(
  object = seurat_merged, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "harmony", dims = 1:30)
seurat_integrated <- FindClusters(seurat_integrated, resolution = 1, cluster.name = "harmony_clusters")
seurat_integrated <- RunUMAP(seurat_integrated, reduction = "harmony", dims = 1:30)

DimPlot(seurat_integrated, label.box = T, label = T)
DimPlot(seurat_integrated, split.by = "condition_time")
#FeaturePlot(seurat_integrated, features = "Tnf", split.by = "condition_time")