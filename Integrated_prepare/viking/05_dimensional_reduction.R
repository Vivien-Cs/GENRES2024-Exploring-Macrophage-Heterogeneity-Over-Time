

#seurat_merged <- FindVariableFeatures(seurat_merged, nfeatures = 1000)
seurat_merged <- SCTransform(seurat_merged)

seurat_merged <- RunPCA(seurat_merged)
#seurat_merged <- RunPCA(seurat_merged, features = gene_names)
ElbowPlot(seurat_merged)
#Normalize & integrate
# options(future.globals.maxSize = 3e+09)
# seurat_merged <- SCTransform(seurat_merged, vars.to.regress = c("percent.mt", "percent.rps", "percent.rpl", "nCount_RNA", "nFeature_RNA"), variable.features.n = 5000)
# seurat_merged <- RunPCA(seurat_merged, npcs = 30, verbose = F)
# seurat_merged <- IntegrateLayers(
#   object = seurat_merged,
#   method = RPCAIntegration,
#   normalization.method = "SCT",
#   verbose = F
# )

# for genes >200 use this
seurat_merged <- FindNeighbors(seurat_merged, reduction = "pca", dims = 1:5)
seurat_merged <- FindClusters(seurat_merged, verbose = FALSE, resolution=res)
seurat_merged <- RunUMAP(seurat_merged, reduction = "pca", dims = 1:5)
seurat_merged <- RunTSNE(seurat_merged, reduction = "pca", dims = 1:5, check_duplicates = F)

# saveRDS(seurat_merged, paste0("preprocessed_merged_dims4_",res,".Rds"))

# looking at clustering using integrated RPC genes <200
# seurat_merged <- FindNeighbors(seurat_merged,features = gene_names)
# seurat_merged <- FindClusters(seurat_merged, resolution = 0.6)
# seurat_merged <- RunUMAP(seurat_merged, features = gene_names)
# seurat_merged <- RunTSNE(seurat_merged, features = gene_names, check_duplicates = F)

DimPlot(seurat_merged, label = TRUE, label.box = TRUE)+NoLegend()+
  DimPlot(seurat_merged, reduction="umap", group.by = "orig.ident",pt.size = 1.5, label.size = 8)+NoLegend()
DimPlot(seurat_merged, reduction="tsne", label = TRUE, label.box = TRUE, pt.size = 1.5, label.size = 8)+NoLegend()
DimPlot(seurat_merged, reduction="tsne", split.by = "time_pt", group.by = "time_pt", cols = c("lightsalmon","magenta4","azure3", "palegreen3","royalblue2"), pt.size = 2)
DimPlot(seurat_merged, reduction="tsne", group.by = "condition", cols = c("cadetblue","brown","goldenrod2"), pt.size = 2)

for(pc in 1:5){
  print(DimPlot(seurat_merged, reduction = "pca", dims = c(pc, pc+1), group.by="time_pt", cols = c("lightsalmon","magenta4","azure3","palegreen3","royalblue2")))
  print(DimPlot(seurat_merged, reduction = "pca", dims = c(pc, pc+1), group.by="condition", cols = c("cadetblue","brown","goldenrod2")))
}