
# for genes >200 use this
seurat_merged <- FindNeighbors(seurat_merged, reduction = "pca", dims = 1:18)
seurat_merged <- FindClusters(seurat_merged, verbose = FALSE, resolution=res)
seurat_merged <- RunUMAP(seurat_merged, reduction = "pca", dims = 1:18)
seurat_merged <- RunTSNE(seurat_merged, reduction = "pca", dims = 1:18, check_duplicates = F)

saveRDS(seurat_merged, paste0("preprocessed_merged_dims4_",res,".Rds"))

#combine time and condition 
seurat_merged$condition_time <- paste0(seurat_merged$condition,"_",seurat_merged$time_pt)

DimPlot(seurat_merged, label = TRUE, label.box = TRUE)+NoLegend()+
  DimPlot(seurat_merged, reduction="umap", group.by = "condition_time",pt.size = 1.5, label.size = 8)+NoLegend()
DimPlot(seurat_merged, reduction="tsne", label = TRUE, label.box = TRUE, pt.size = 1.5, label.size = 8)+NoLegend()
DimPlot(seurat_merged, reduction="tsne", split.by = "time_pt", group.by = "time_pt", cols = c("lightsalmon","magenta4","azure3", "palegreen3","royalblue2"), pt.size = 2)
DimPlot(seurat_merged, reduction="tsne", group.by = "condition", cols = c("cadetblue","brown","goldenrod2"), pt.size = 2)

for(pc in 1:5){
  print(DimPlot(seurat_merged, reduction = "pca", dims = c(pc, pc+1), group.by="time_pt", cols = c("lightsalmon","magenta4","azure3","palegreen3","royalblue2")))
  print(DimPlot(seurat_merged, reduction = "pca", dims = c(pc, pc+1), group.by="condition", cols = c("cadetblue","brown","goldenrod2")))
}
saveRDS(pca_plots, "pca_plots.rds")