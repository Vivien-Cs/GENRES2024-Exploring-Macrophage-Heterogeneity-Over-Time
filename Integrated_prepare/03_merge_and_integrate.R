### merge datasets ###
seurat_merged <- merge(
  seuratObjList[[1]], 
  y = seuratObjList[2:length(seuratObjList)]
)

#merge by condition and time 
seurat_merged$condition_time <- paste0(as.character(seurat_merged$condition), "_", as.character(seurat_merged$time_pt))

# Reorder the levels of the condition_time factor to match the desired order
seurat_merged$condition_time <- factor(seurat_merged$condition_time, 
                                       levels = c("M0_Liu_2020_zero","Media_48_Li_2019_zero","Media_M0_Hagai_2018_zero", "M0_Liu_2020_six", "M0_Liu_2020_twenty4", 
                                                  "LPS_M1_Hagai_2018_two", "LPS_M1_Hagai_2018_four", "LPS_M1_Hagai_2018_six", 
                                                  "M1_LPS_IFNG_Liu_2020_six","M1_LPS_IFNG_Liu_2020_twenty4",
                                                  "M1_LPS_IFNG_48_Li_2019_fourty8"))

#SCTransform setting method as glmGamPoi
seurat_merged <- SCTransform(seurat_merged, vst.flavor = "v2", do.center = T, vars.to.regress = c("percent.mt","percent.rps","percent.rpl"))

top_2000 <- head(VariableFeatures(seurat_merged), 2000)
write.csv(top_2000, paste0(output_integrated,"top_2000.csv"), row.names = FALSE)

#Run PCA - Standard Workflow
seurat_merged <- ScaleData(seurat_merged, features = top_2000)
seurat_merged_standard <- RunPCA(seurat_merged, features = top_2000)

ElbowPlot(seurat_merged_standard) + ggtitle(label = "Standard Workflow")

#Run integration#
seurat_integrated_standard <- IntegrateLayers(
  object = seurat_merged_standard, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  normalization.method = "SCT",
  verbose = FALSE 
)

# Save integrated object
saveRDS(seurat_integrated, file = paste0(output_integrated, "seurat_integrated.rds"))  

##################################
##################################
##################################

########### Run PCA - AMDSG Workflow ###########

# Load data based on the AMDSG and PSG dataset
AMDSG <- read_excel("G:/My Drive/Viv_MACHS/R/AMDSG_geneID.xlsx") 
amdsg_gene <- AMDSG$Gene
PSG <- read_excel("G:/My Drive/Viv_MACHS/R/PSG_geneID.xlsx")
psg_gene <- PSG$Gene
gene_names <- unique(c(amdsg_gene, psg_gene)
)

# Scale all genes in gene_names
seurat_merged_genelist <- ScaleData(seurat_merged, features = gene_names, verbose = FALSE)

seurat_merged_genelist <- RunPCA(seurat_merged, features = gene_names, verbose = FALSE)

#merge by condition and time
seurat_merged_genelist$condition_time <- paste0(as.character(seurat_merged_standard$condition), "_", as.character(seurat_merged_standard$time_pt))

# Reorder the levels of the condition_time factor to match the desired order
seurat_merged_genelist$condition_time <- factor(seurat_merged_genelist$condition_time, 
                                       levels = c("M0_Liu_2020_zero","Media_48_Li_2019_zero","Media_M0_Hagai_2018_zero", "M0_Liu_2020_six", "M0_Liu_2020_twenty4", 
                                                  "LPS_M1_Hagai_2018_two", "LPS_M1_Hagai_2018_four", "LPS_M1_Hagai_2018_six", 
                                                  "M1_LPS_IFNG_Liu_2020_six","M1_LPS_IFNG_Liu_2020_twenty4",
                                                  "M1_LPS_IFNG_48_Li_2019_fourty8"))

ElbowPlot(seurat_merged_genelist) + ggtitle(label = "AMDSG_PSG Workflow")

#Run integration#
seurat_integrated_genelist <- IntegrateLayers(
  object = seurat_merged_genelist, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  normalization.method = "SCT",
  verbose = FALSE 
)

# Save integrated object
saveRDS(seurat_integrated_genelist, file = paste0(output_integrated, "seurat_integrated_amdsg_psg.rds"))
##################################
##################################
##################################


rm(seurat_merged_standard, seurat_merged_genelist, seurat_merged)


