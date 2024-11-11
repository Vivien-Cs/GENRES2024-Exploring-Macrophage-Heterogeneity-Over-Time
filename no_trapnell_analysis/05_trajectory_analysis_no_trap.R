
#seurat_integrated_genelist – This is based on the AMDSG workflow.
#seurat_integrated_standard – This is based on the standard workflow.

# Extract PCA embeddings and clusters from the integrated object
rd <- Embeddings(seurat_integrated_genelist, "pca")[, 1:5] 
cl <- seurat_integrated_genelist$seurat_clusters

# Convert to matrix and factor as required by Slingshot
rd <- as.matrix(rd)
cl <- as.factor(cl)

# Run Slingshot to infer trajectories
set.seed(40)
slingshot_obj <- slingshot(rd, clusterLabels = cl)

# Identify cells from the starting samples
starting_samples <- c("Med_M0_0h_rep1", "Med_M0_0h_rep2", "Med_M0_0h_rep3", "M0_0h_48", "M0_0h")

# Extract the metadata from the Seurat object to verify sample assignments
cell_metadata <- seurat_integrated_genelist@meta.data

# Check the cells that belong to the start clusters/samples
starting_cells <- which(seurat_integrated_genelist$orig.ident %in% starting_samples)

# Use the clusters from the starting cells to determine the starting clusters
start_cluster <- unique(cl[starting_cells])

# Running Slingshot with the start clusters
slingshot_obj <- slingshot(rd, clusterLabels = cl, start.clus = start_cluster)


                                           #### pseudotime ##
# Extract pseudotime from the Slingshot object
pseudotime <- slingPseudotime(slingshot_obj)

# Add pseudotime as metadata to the Seurat object
seurat_integrated_genelist$pseudotime <- pseudotime

# Visualise UMAP with pseudotime
p1 <- DimPlot(seurat_integrated_genelist, reduction = "umap", group.by = "pseudotime", label = TRUE) +
  scale_color_viridis_c(option = "plasma") +
  labs(title = "UMAP with Pseudotime", x = "UMAP1", y = "UMAP2") +
  theme_minimal()
# Save the plot
ggsave("G:/My Drive/Viv_MACHS/R/markers_merged/UMAP_with_Pseudotime_slingshot.png", plot = p1, width = 8, height = 6, dpi = 300)





                           ########### monocle3 ##############
# 2. Convert Seurat object to Monocle 3 CDS
expression_matrix <- GetAssayData(seurat_integrated_genelist, layer = "data")

cell_metadata <- seurat_integrated_genelist@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix), row.names = rownames(expression_matrix))

# Creating the Monocle 3 CDS object manually
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

# 3. Preprocess the CDS and reduce dimensions using PCA or UMAP
cds <- preprocess_cds(cds, method = "PCA", num_dim = 30)  # Ensure PCA is set
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA", umap.n_neighbors = 30, umap.min_dist = 0.1)

# 4. Cluster cells 
cds <- cluster_cells(cds)

# 5. Learn the trajectory graph
cds <- learn_graph(cds)

# 6. Identify starting cells based on your root samples
starting_samples <- c("Med_M0_0h_rep1", "Med_M0_0h_rep2", "Med_M0_0h_rep3", "M0_0h_48", "M0_0h")
root_cells <- colnames(cds)[cds@colData$orig.ident %in% starting_samples]

# 7. Order cells in pseudotime starting from the specified root cells
cds <- order_cells(cds, root_cells = root_cells)

# 8. Check pseudotime values to ensure they are calculated correctly
pseudotime_values <- pseudotime(cds)
head(pseudotime_values)

# 9. Identify genes significantly associated with pseudotime
cds_pr_test_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
significant_genes <- subset(cds_pr_test_res, q_value < 0.05)$gene_short_name

# 10. Visualize pseudotime in UMAP
p2 <- plot_cells(cds,
           color_cells_by = "pseudotime",      # Color by pseudotime
           show_trajectory_graph = TRUE,       # Show the trajectory graph
           trajectory_graph_color = "red",     # Make the trajectory graph more visible
           label_groups_by_cluster = TRUE,     # Label clusters
           label_leaves = TRUE,                # Label leaves
           label_branch_points = TRUE,         # Label branch points
           graph_label_size = 5,               # Increase size of labels
           cell_size = 1.8)                    # Adjust cell size for better visibility
# Save the UMAP pseudotime trajectory plot
ggsave("G:/My Drive/Viv_MACHS/R/markers_merged/UMAP_with_trajectory.png", plot = p2, width = 8, height = 6, dpi = 300)

# 11. Plot heatmap of significant genes along pseudotime
p3 <- plot_pseudotime_heatmap(cds[significant_genes, ], 
                        num_clusters = 4,  # Set the number of gene clusters
                        show_rownames = TRUE)  # Show gene names on y-axis

# Save the heatmap
ggsave("G:/My Drive/Viv_MACHS/R/markers_merged/Pseudotime_Heatmap.png", plot = p3, width = 8, height = 6, dpi = 300)