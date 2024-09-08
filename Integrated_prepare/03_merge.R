### merge datasets ###
seurat_merged <- merge(
  seuratObjList[[lstSamples[[1]]]], 
  y = seuratObjList[2:length(seuratObjList)]
)

#SCTransform setting method as glmGamPoi
seurat_merged <- SCTransform(seurat_merged, vst.flavor = "v2")

#top_2000 <- head(VariableFeatures(seurat_merged), 2000)
#write.csv(top_2000, "G:/My Drive/Genres2024/output/top_2000.csv", row.names = FALSE)

#COME BACK TO THIS!!!
# variable_plot <-VariableFeatSurePlot(
#   seurat_merged,
#   cols = c("black", "red"),
#   pt.size = 1,
# ) + LabelPoints(plot = VariableFeaturePlot(seurat_merged),
#                 points = top_500,
#                 repel = TRUE)
# print(variable_plot)

seurat_merged <- ScaleData(seurat_merged, features = gene_names)

seurat_merged <- RunPCA(seurat_merged, features = gene_names)

ElbowPlot(seurat_merged)


# # SCTransform and merging can be saved as an RDS file for faster run, outline a path for it to save
# output_path <- "G:/My Drive/Genres2024/R/"

# # Save the Seurat object using the new path
# saveRDS(seurat_merged, file = paste0(output_path, "processed_seurat_merged_dims", dims, "_res", res, ".rds"))

#merge by condition and time
seurat_merged$condition_time <- paste0(seurat_merged$condition,"_",seurat_merged$time_pt)

# Reorder the levels of the condition_time factor to match the desired order
seurat_merged$condition_time <- factor(seurat_merged$condition_time, 
                                      levels = c("M0_Liu_2020_zero","Media_48_JCI_zero","Media_M0_Hagai_2018_zero", "M0_Liu_2020_six", "M0_Liu_2020_twenty4", 
                                                 "LPS_M1_Hagai_2018_two", "LPS_M1_Hagai_2018_four", "LPS_M1_Hagai_2018_six", 
                                                 "M1_LPS_IFNG_Liu_2020_six","M1_LPS_IFNG_Liu_2020_twenty4",
                                                 "M1_LPS_IFNG_48_JCI_fourty8"))
DimPlot(seurat_merged, group.by = "condition")
# Save the plot as a PDF file
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

# seurat_merged$paper <- factor (seurat_merged$paper, levels = c("Liu", "Hag", "JCI"))
# Dimplo(Seurat_merged, group.by = "paper", label = TRUE) +
#   ggtitle("UMAP grouped by paper")