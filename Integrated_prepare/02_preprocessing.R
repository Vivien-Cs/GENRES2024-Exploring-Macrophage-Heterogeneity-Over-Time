
                      ########   perform standard preprocessing on each object   ###########
#check gene types
#grep(pattern = "^Ccl", x = rownames(seuratObj), value = TRUE)
for (i in 1:length(seuratObjList)) {
  #seuratObjList[[i]]
  #'Gm42418 ','AY036118' removed from analysis for finding markers
  #"Genes Gm42418 and AY036118 were also removed, as they overlap the rRNA element Rn45s and represent rRNA contamination"
  #Source:https://www.nature.com/articles/s41467-021-27035-8
  counts <- GetAssayData(seuratObjList[[i]], assay = "RNA")
  counts <- counts[-(which(rownames(counts) %in% c('Gm42418','AY036118'))),]
  seuratObjList[[i]] <- subset(seuratObjList[[i]], features = rownames(counts))
  seuratObjList[[i]] <- PercentageFeatureSet(seuratObjList[[i]], pattern = "^mt-", col.name = "percent.mt")
  seuratObjList[[i]] <- subset(seuratObjList[[i]], subset = percent.mt < 10)
  seuratObjList[[i]] <- PercentageFeatureSet(seuratObjList[[i]], pattern = "^Rps", col.name = "percent.rps")
  seuratObjList[[i]] <- PercentageFeatureSet(seuratObjList[[i]], pattern = "^Rpl", col.name = "percent.rpl")
  print(seuratObjList[[i]])
  saveRDS(seuratObjList[[i]], file = paste0(output_integrated, "preprocessed_seuratObj_", lstSamples[[i]], ".rds")) 
}

# seuratObjList <- lapply(X = seuratObjList, FUN = SCTransform, method = "glmGamPoi", vars.to.regress = c("percent.mt", "percent.rps", "percent.rpl", "nCount_RNA", "nFeature_RNA"), verbose = TRUE)