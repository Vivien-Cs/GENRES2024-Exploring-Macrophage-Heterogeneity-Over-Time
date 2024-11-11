########### TOP GENES IN CLUSTER 1 AND 3 ###########
# Identify Top Marker Genes for Clusters 1 and 3

#finding markers with SCTransform
seurat_integrated_genelist <- PrepSCTFindMarkers(seurat_integrated_genelist)

# Find markers for clusters 1 and 3
cluster1_genes <- FindMarkers(seurat_integrated_genelist, ident.1 = 1, logfc.threshold = 0.25)
cluster3_genes <- FindMarkers(seurat_integrated_genelist, ident.1 = 3, logfc.threshold = 0.25)

# Select top 50 genes for each cluster
top_genes_cluster1 <- rownames(cluster1_genes)[1:100]
top_genes_cluster3 <- rownames(cluster3_genes)[1:100]

# Extract and prepare expression data over time

# Function to prepare data for line plots
prepare_line_plot_data <- function(seurat_integrated_genelist, top_genes, cluster_label) {
  gene_expr_data <- FetchData(seurat_integrated_genelist, vars = c("time_pt", top_genes)) %>%
    group_by(time_pt) %>%
    summarize(across(all_of(top_genes), mean)) %>%
    pivot_longer(cols = -time_pt, names_to = "gene", values_to = "expression") %>%
    mutate(cluster = cluster_label)
  
  return(gene_expr_data)
}

# Prepare data for cluster 1
gene_expr_long_cluster1 <- prepare_line_plot_data(seurat_integrated_genelist, top_genes_cluster1, "Cluster 1")

# Prepare data for cluster 3
gene_expr_long_cluster3 <- prepare_line_plot_data(seurat_integrated_genelist, top_genes_cluster3, "Cluster 3")

# Plot Gene Expression Over Time for Each Cluster Separately

# Function to create line plots
plot_gene_expression_over_time <- function(gene_expr_long_data, cluster_label) {
  ggplot(gene_expr_long_data, aes(x = time_pt, y = expression, color = gene, group = gene)) +
    geom_line() +
    facet_wrap(~ gene, scales = "free_y") +
    labs(title = paste("Gene Expression Over Time for", cluster_label), x = "Time Point", y = "Expression") +
    theme_minimal() +
    theme(legend.position = "none")
}

# Plot for cluster 1
plot_cluster1 <- plot_gene_expression_over_time(gene_expr_long_cluster1, "Cluster 1")
ggsave(filename = "G:/My Drive/Viv_MACHS/R/plot_cluster1_top_100_sep.png", plot = plot_cluster1, width = 10, height = 8, dpi = 300 , bg = "white")

# Single graph plot for gene expression over time with all genes on the same axes
cluster_1_top_100 <- ggplot(gene_expr_long_cluster1, aes(x = time_pt, y = expression, color = gene, group = gene)) +
  geom_line(size = 1) +               # Add lines for each gene
  geom_point(size = 2) +              # Optionally add points for each time point
  labs(title = "Gene Expression Over Time for Cluster 1", x = "Time Point", y = "Expression") +
  scale_color_discrete(name = "Gene") +  # Label the legend as "Gene"
  theme_minimal() +
  theme(legend.position = "right")

ggsave(filename = "G:/My Drive/Viv_MACHS/R/cluster_1_top_100.png", plot = cluster_1_top_100, width = 10, height = 8, dpi = 300, bg = "white")

# Plot for cluster 3
plot_cluster3 <- plot_gene_expression_over_time(gene_expr_long_cluster3, "Cluster 3")
print(plot_cluster3)
ggsave(filename = "G:/My Drive/Viv_MACHS/R/plot_cluster3_top_100_sep.png", plot = plot_cluster3, width = 10, height = 8, dpi = 300, bg = "white")

# single graph plot for cluster 3
cluster_3_top_100 <- ggplot(gene_expr_long_cluster3, aes(x = time_pt, y = expression, color = gene, group = gene)) +
  geom_line(size = 1) +             
  geom_point(size = 2) +           
  labs(title = "Gene Expression Over Time for Cluster 3", x = "Time Point", y = "Expression") +
  scale_color_discrete(name = "Gene") +  
  theme_minimal() +
  theme(legend.position = "right")     

ggsave(filename = "G:/My Drive/Viv_MACHS/R/cluster_3_top_100.png", plot = cluster_3_top_100, width = 10, height = 8, dpi = 300, bg = "white")

####################### addmodule score ######################################
# Define the circuit with the specified genes
ccl_circuit_genes <- list(CCL_Circuit = c("Ccl5", "Ccl8", "Ccl1"))

# Add module score to the Seurat object for the CCL circuit
seurat_integrated_genelist <- AddModuleScore(seurat_integrated_genelist, features = ccl_circuit_genes, name = "CCL_CircuitScore")

# The module score will be added as a column named "CCL_CircuitScore1" in the Seurat object metadata
# Violin plot to see distribution of the module score across clusters
VlnPlot(seurat_integrated_genelist, features = "CCL_CircuitScore1", group.by = "seurat_clusters") +
  labs(title = "CCL Circuit Score Across Clusters")

# Feature plot to see spatial distribution on UMAP or other embedding
FeaturePlot(seurat_integrated_genelist, features = "CCL_CircuitScore1") +
  labs(title = "CCL Circuit Score on UMAP")

# Fetch data for the time points and the CCL circuit score
ccl_module_scores_data <- FetchData(seurat_integrated_genelist, vars = c("time_pt", "CCL_CircuitScore1"))

# Calculate mean module score for each time point
ccl_module_scores_avg <- ccl_module_scores_data %>%
  group_by(time_pt) %>%
  summarize(CCL_CircuitScore = mean(CCL_CircuitScore1, na.rm = TRUE))

# Plot module score over time points
ggplot(ccl_module_scores_avg, aes(x = time_pt, y = CCL_CircuitScore)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue", size = 2) +
  labs(title = "CCL Circuit Score Over Time", x = "Time Point", y = "CCL Circuit Score") +
  theme_minimal()



####### add module for cluster 1 and 3 #########
#Add the module score (if not already done)
ccl_circuit_genes <- list(CCL_Circuit = c("Ccl5", "Ccl8", "Ccl1"))
seurat_integrated_genelist <- AddModuleScore(seurat_integrated_genelist, features = ccl_circuit_genes, name = "CCL_CircuitScore")

#Fetch the data for CCL circuit score, cluster, and time point
ccl_module_scores_data <- FetchData(seurat_integrated_genelist, vars = c("time_pt", "seurat_clusters", "CCL_CircuitScore1"))

# Filter the data for cluster 1 & cluster 3 and calculate mean scores over time

# For Cluster 1
ccl_module_scores_cluster1 <- ccl_module_scores_data %>%
  filter(seurat_clusters == 1) %>%
  group_by(time_pt) %>%
  summarize(CCL_CircuitScore = mean(CCL_CircuitScore1, na.rm = TRUE))

# For Cluster 3
ccl_module_scores_cluster3 <- ccl_module_scores_data %>%
  filter(seurat_clusters == 3) %>%
  group_by(time_pt) %>%
  summarize(CCL_CircuitScore = mean(CCL_CircuitScore1, na.rm = TRUE))

#Plot the CCL circuit score over time for each cluster separately

# Plot for Cluster 1
plot_cluster1 <- ggplot(ccl_module_scores_cluster1, aes(x = time_pt, y = CCL_CircuitScore)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue", size = 2) +
  labs(title = "CCL Circuit Score Over Time in Cluster 1", x = "Time Point", y = "CCL Circuit Score") +
  theme_minimal()

# Plot for Cluster 3
plot_cluster3 <- ggplot(ccl_module_scores_cluster3, aes(x = time_pt, y = CCL_CircuitScore)) +
  geom_line(color = "red", size = 1) +
  geom_point(color = "red", size = 2) +
  labs(title = "CCL Circuit Score Over Time in Cluster 3", x = "Time Point", y = "CCL Circuit Score") +
  theme_minimal()

# Display the plots
print(plot_cluster1)
print(plot_cluster3)

# Save 
ggsave(filename = "G:/My Drive/Viv_MACHS/R/CCL_Circuit_Cluster1.png", plot = plot_cluster1, width = 8, height = 6, dpi = 300, bg = "white")
ggsave(filename = "G:/My Drive/Viv_MACHS/R/CCL_Circuit_Cluster3.png", plot = plot_cluster3, width = 8, height = 6, dpi = 300, bg = "white")




################ check ccl individually ########
#get expression data for each gene in the circuit
gene_data <- FetchData(seurat_integrated_genelist, vars = c("time_pt", "seurat_clusters", "Ccl5", "Ccl8", "Ccl1"))

#  Filter for Cluster 1 and Cluster 3, then calculate mean expression for each gene over time

# For Cluster 1
gene_data_cluster1 <- gene_data %>%
  filter(seurat_clusters == 1) %>%
  group_by(time_pt) %>%
  summarize(Ccl5 = mean(Ccl5, na.rm = TRUE), 
            Ccl8 = mean(Ccl8, na.rm = TRUE), 
            Ccl1 = mean(Ccl1, na.rm = TRUE))

# For Cluster 3
gene_data_cluster3 <- gene_data %>%
  filter(seurat_clusters == 3) %>%
  group_by(time_pt) %>%
  summarize(Ccl5 = mean(Ccl5, na.rm = TRUE), 
            Ccl8 = mean(Ccl8, na.rm = TRUE), 
            Ccl1 = mean(Ccl1, na.rm = TRUE))

# Plot each gene separately over time in Cluster 1 and Cluster 3

# For Cluster 1
plot_cluster1_genes <- gene_data_cluster1 %>%
  pivot_longer(cols = c("Ccl5", "Ccl8", "Ccl1"), names_to = "gene", values_to = "expression") %>%
  ggplot(aes(x = time_pt, y = expression, color = gene, group = gene)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Expression of CCL Circuit Genes Over Time in Cluster 1", x = "Time Point", y = "Expression") +
  theme_minimal()

# For Cluster 3
plot_cluster3_genes <- gene_data_cluster3 %>%
  pivot_longer(cols = c("Ccl5", "Ccl8", "Ccl1"), names_to = "gene", values_to = "expression") %>%
  ggplot(aes(x = time_pt, y = expression, color = gene, group = gene)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Expression of CCL Circuit Genes Over Time in Cluster 3", x = "Time Point", y = "Expression") +
  theme_minimal()

# Display the plots
print(plot_cluster1_genes)
print(plot_cluster3_genes)

# save
ggsave(filename = "G:/My Drive/Viv_MACHS/R/CCL_Genes_Cluster1.png", plot = plot_cluster1_genes, width = 8, height = 6, dpi = 300, bg = "white")
ggsave(filename = "G:/My Drive/Viv_MACHS/R/CCL_Genes_Cluster3.png", plot = plot_cluster3_genes, width = 8, height = 6, dpi = 300, bg = "white")



################# Check for TNF & NOS2 #########
# Get expression data for TNF and NOS2 along with cluster and time point information
tnf_nos2_data <- FetchData(seurat_integrated_genelist, vars = c("time_pt", "seurat_clusters", "Tnf", "Nos2"))

# Filter for Cluster 1 and Cluster 3, then calculate mean expression for TNF and NOS2 over time

# For Cluster 1
tnf_nos2_data_cluster1 <- tnf_nos2_data %>%
  filter(seurat_clusters == 1) %>%
  group_by(time_pt) %>%
  summarize(Tnf = mean(Tnf, na.rm = TRUE), 
            Nos2 = mean(Nos2, na.rm = TRUE))

# For Cluster 3
tnf_nos2_data_cluster3 <- tnf_nos2_data %>%
  filter(seurat_clusters == 3) %>%
  group_by(time_pt) %>%
  summarize(Tnf = mean(Tnf, na.rm = TRUE), 
            Nos2 = mean(Nos2, na.rm = TRUE))

# Plot TNF and NOS2 separately over time in Cluster 1 and Cluster 3

# For Cluster 1
plot_tnf_nos2_cluster1 <- tnf_nos2_data_cluster1 %>%
  pivot_longer(cols = c("Tnf", "Nos2"), names_to = "gene", values_to = "expression") %>%
  ggplot(aes(x = time_pt, y = expression, color = gene, group = gene)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Expression of TNF and NOS2 Over Time in Cluster 1", x = "Time Point", y = "Expression") +
  theme_minimal()

# For Cluster 3
plot_tnf_nos2_cluster3 <- tnf_nos2_data_cluster3 %>%
  pivot_longer(cols = c("Tnf", "Nos2"), names_to = "gene", values_to = "expression") %>%
  ggplot(aes(x = time_pt, y = expression, color = gene, group = gene)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Expression of TNF and NOS2 Over Time in Cluster 3", x = "Time Point", y = "Expression") +
  theme_minimal()

# Display the plots
print(plot_tnf_nos2_cluster1)
print(plot_tnf_nos2_cluster3)

# Save
ggsave(filename = "G:/My Drive/Viv_MACHS/R/TNF_NOS2_Cluster1.png", plot = plot_tnf_nos2_cluster1, width = 8, height = 6, dpi = 300, bg = "white")
ggsave(filename = "G:/My Drive/Viv_MACHS/R/TNF_NOS2_Cluster3.png", plot = plot_tnf_nos2_cluster3, width = 8, height = 6, dpi = 300, bg = "white")


