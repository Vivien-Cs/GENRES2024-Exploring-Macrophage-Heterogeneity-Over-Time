##integrate using entire dataset##
# seurat_merged <- readRDS(file = paste0(output_path, "processed_seurat_merged_dims", dims, "_res", res, ".rds"))

seurat_integrated <- IntegrateLayers(
  object = seurat_merged, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "harmony", dims = 1:10)
seurat_integrated <- FindClusters(seurat_integrated, resolution = 0.2, cluster.name = "harmony_clusters")
seurat_integrated <- RunUMAP(seurat_integrated, reduction = "harmony", dims = 1:10)

all_cell <- DimPlot(seurat_integrated, label.box = T, label = T)

# Save the plot as a PDF file
ggsave("G:/My Drive/Genres2024/output/all_cells_0.2_standard.pdf", plot = all_cell, width = 10, height = 8)

# Reorder the levels of the time_pt factor
seurat_integrated$time_pt <- factor(seurat_integrated$time_pt, levels = c("zero", "two", "four", "six", "twenty4", "fourty8"))


# Create the plot and store it in a variable
split_plot <- DimPlot(seurat_integrated, split.by = "condition", ncol= 5)
# Save the plot as a PDF file
ggsave("G:/My Drive/Genres2024/output/split_by_condition_0.2_standard.pdf", plot = split_plot, width = 35, height = 10)

#run tnf
tnf_split_plot <- FeaturePlot(seurat_integrated, features = "Tnf", split.by = "condition_time", ncol = 5)
ggsave("G:/My Drive/Genres2024/output/tnf_by_condition_time_standard_0.02.pdf", plot = tnf_split_plot, width = 30, height = 20)



######### Count table ################
condition_counts <- table(Idents(seurat_integrated), seurat_integrated$condition)
print(condition_counts)

# Percentage table
condition_percentages <- prop.table(condition_counts, margin = 2) * 100
print(round(condition_percentages, 2))

# Save to CSV
write.csv(condition_counts, "G:/My Drive/Genres2024/output/condition_counts_standard.csv")
write.csv(condition_percentages, "G:/My Drive/Genres2024/output/condition_percentages_standard.csv")

## condition_Time
# Count table
condition_time_counts <- table(Idents(seurat_integrated), seurat_integrated$condition_time)
print(condition_time_counts)

# Percentage table
condition_time_percentages <- prop.table(condition_time_counts, margin = 2) * 100
print(round(condition_time_percentages, 2))

# Save to CSV
write.csv(condition_time_counts, "G:/My Drive/Genres2024/output/condition_time_counts_standard.csv")
write.csv(condition_time_percentages, "G:/My Drive/Genres2024/output/condition_time_percentages_standard.csv")

# Count table
time_counts <- table(Idents(seurat_integrated), seurat_integrated$time_pt)
print(time_counts)

# Percentage table
time_percentages <- prop.table(time_counts, margin = 2) * 100
print(round(time_percentages, 2))

# Save to CSV
write.csv(time_counts, "G:/My Drive/Genres2024/output/time_counts_standard.csv")
write.csv(time_percentages, "G:/My Drive/Genres2024/output/time_percentages_standard.csv")

#FeaturePlot(seurat_integrated, features = "Tnf", split.by = "condition_time")