### merge datasets ###
seurat_merged <- merge(
  seuratObjList[[lstSamples[[1]]]], 
  y = seuratObjList[2:length(seuratObjList)]
)
#SCTransform setting method as glmGamPoi
seurat_merged <- SCTransform(seurat_merged, vst.flavor = "v2")

seurat_merged <- RunPCA(seurat_merged)
ElbowPlot(seurat_merged)

saveRDS(seurat_merged, file = paste0(output_integrated, "processed_seurat_merged_dims", dims, "_res", res, ".rds"))

#merge by condition and time
seurat_merged$condition_time <- paste0(seurat_merged$condition,"_",seurat_merged$time_pt)