
                      ########   perform standard preprocessing on each object   ###########
#check gene types
#grep(pattern = "^Ccl", x = rownames(seuratObj), value = TRUE)
all_samples <- c(lstSamples, lstSamples_48h, lstSamples_2_4_6h)
for (sample_name in names(seuratObjList)) {
  seuratObj <- seuratObjList[[sample_name]]
  
  # Perform preprocessing steps
  counts <- GetAssayData(seuratObj, assay = "RNA")
  counts <- counts[!rownames(counts) %in% c('Gm42418','AY036118'),]
  seuratObj <- subset(seuratObj, features = rownames(counts))
  seuratObj <- PercentageFeatureSet(seuratObj, pattern = "^mt-", col.name = "percent.mt")
  seuratObj <- subset(seuratObj, subset = percent.mt < 10)
  seuratObj <- PercentageFeatureSet(seuratObj, pattern = "^Rps", col.name = "percent.rps")
  seuratObj <- PercentageFeatureSet(seuratObj, pattern = "^Rpl", col.name = "percent.rpl")
  
  print(seuratObj)
  
  # Update the object in the list
  seuratObjList[[sample_name]] <- seuratObj
  
  # Save the preprocessed object
  saveRDS(seuratObj, file = paste0(output_integrated, "preprocessed_seuratObj_", sample_name, ".rds"))
}
# seuratObjList <- lapply(X = seuratObjList, FUN = SCTransform, method = "glmGamPoi", vars.to.regress = c("percent.mt", "percent.rps", "percent.rpl", "nCount_RNA", "nFeature_RNA"), verbose = TRUE)
                      
                      