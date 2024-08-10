### merge datasets
seurat_merged <- merge(
  seuratObjList[[lstSamples[[1]]]], 
  y = seuratObjList[2:length(seuratObjList)]
)


