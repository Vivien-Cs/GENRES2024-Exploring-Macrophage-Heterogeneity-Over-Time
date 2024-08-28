#Run AMDSG list integration#
seurat_integrated <- IntegrateLayers(
  object = seurat_merged, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

valid_genes <- gene_names[gene_names %in% rownames(seurat_integrated)]

#seurat_integrated <- RunUMAP(seurat_integrated, features = valid_genes)
ElbowPlot(seurat_integrated)

seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:15, reduction = "harmony")
seurat_integrated <- FindClusters(seurat_integrated, resolution = 0.2, cluster.name = "AMDSG_48hr_clusters")
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:15, reduction = "harmony")

#RUN ALL CELLS
all_cell <- DimPlot(seurat_integrated, label.box = T, label = T)
#SAVE
ggsave("G:/My Drive/Genres2024/output/all_cells_0.02_amdsg_psg.pdf", plot = all_cell, width = 10, height = 8)

# Reorder the levels of the time_pt factor
seurat_integrated$time_pt <- factor(seurat_integrated$time_pt, levels = c("zero", "two", "four", "six", "twenty4", "fourty8"))
# seurat_integrated$condition <- factor(seurat_integrated$condition, 
#                                       levels = c("M0_Liu_2020", "M1_LPS_IFNG_Liu_2020", "Media_48_JCI", "M1_LPS_IFNG_48_JCI", "Media_M0_Hagai_2018", "LPS_M1_Hagai_2018"))

split_plot <- DimPlot(seurat_integrated, split.by = "condition_time", ncol = 5)
ggsave("G:/My Drive/Genres2024/output/split_by_condition_time_amdsg_psg_0.02.pdf", plot = split_plot, width = 15, height = 20)

tnf_split_plot <- FeaturePlot(seurat_integrated, features = "Tnf", split.by = "condition_time", ncol = 4)
ggsave("G:/My Drive/Genres2024/output/tnf_by_condition_time_amdsg_0.02.pdf", plot = tnf_split_plot, width = 30, height = 20)



### count###

# Count table
condition_counts <- table(Idents(seurat_integrated), seurat_integrated$condition)
print(condition_counts)

# Percentage table
condition_percentages <- prop.table(condition_counts, margin = 2) * 100
print(round(condition_percentages, 2))

# Save to CSV
write.csv(condition_counts, "G:/My Drive/Genres2024/output/condition_counts_amdsg.csv")  ##ACCIDENTALLY SAVED AS STANDARD SO RE-RUN STANDARD HERE#
write.csv(condition_percentages, "G:/My Drive/Genres2024/output/condition_percentages_amdsg.csv")

## condition_Time
# Count table
condition_time_counts <- table(Idents(seurat_integrated), seurat_integrated$condition_time)
print(condition_time_counts)

# Percentage table
condition_time_percentages <- prop.table(condition_time_counts, margin = 2) * 100
print(round(condition_time_percentages, 2))

# Save to CSV
write.csv(condition_time_counts, "G:/My Drive/Genres2024/output/condition_time_counts_amdsg.csv")
write.csv(condition_time_percentages, "G:/My Drive/Genres2024/output/condition_time_percentages_amdsg.csv")

# Count table
time_counts <- table(Idents(seurat_integrated), seurat_integrated$time_pt)
print(time_counts)

# Percentage table
time_percentages <- prop.table(time_counts, margin = 2) * 100
print(round(time_percentages, 2))

# Save to CSV
write.csv(time_counts, "G:/My Drive/Genres2024/output/time_counts_amdsg.csv")
write.csv(time_percentages, "G:/My Drive/Genres2024/output/time_percentages_amdsg.csv")

