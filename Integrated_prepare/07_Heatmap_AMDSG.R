# Extract the data from the Seurat object
sct_data <- GetAssayData(seurat_integrated, layer = "data", assay = "SCT")

# Subset the data to include only the genes of interest (AMDSG genes)
gene_data <- sct_data[rownames(sct_data) %in% gene_names, ]

# Scale the data for better visualisation
scaled_data <- t(scale(t(gene_data)))

# Define annotation for the heatmap
annotation <- seurat_integrated@meta.data[, c("AMDSG_48hr_clusters", "condition", "time_pt")]

# Create the heatmap
ht <- Heatmap(
  matrix = scaled_data, 
  name = "Expression",
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  top_annotation = HeatmapAnnotation(df = annotation)
)

# Save the heatmap as a PNG file
ggsave(
  filename = "G:/My Drive/Genres2024/heatmap_results/AMDSG_gene_heatmap.png", 
  plot = grid.grabExpr(draw(ht)), 
  device = "png", 
  width = 10, 
  height = 8, 
  units = "in", 
  dpi = 300
)

print(heatmap_plot)

# or save as pdf
pdf(filename = "G:/My Drive/Genres2024/heatmap_results/AMDSG_gene_heatmap.pdf", 
    width = 10, height = 8)
