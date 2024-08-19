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

# ######## total cell count ###############
# # Create a table of cell counts by condition and time point
# cell_counts <- table(seurat_merged$condition_time)
# 
# # Convert the table to a data frame
# cell_counts_df <- as.data.frame(cell_counts)
# 
# # Display the table
# print(cell_counts_df)
# # Specify the output file path
# output_file <- "G:/My Drive/Genres2024/output/cell_counts_by_condition_time.csv"
# 
# # Save the table to a CSV file
# write.csv(cell_counts_df, file = output_file, row.names = FALSE)
