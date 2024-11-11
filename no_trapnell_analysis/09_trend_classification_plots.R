                          ##### sort by trend CLUSTER 1########
#Identify the top marker genes for Cluster 1
seurat_integrated_genelist <- PrepSCTFindMarkers(seurat_integrated_genelist)
cluster1_genes <- FindMarkers(seurat_integrated_genelist, ident.1 = 1, logfc.threshold = 0.25)
top_genes_cluster1 <- rownames(cluster1_genes)[1:100]

# Prepare data for line plot
gene_expr_data <- FetchData(seurat_integrated_genelist, vars = c("time_pt", top_genes_cluster1)) %>%
  group_by(time_pt) %>%
  summarize(across(all_of(top_genes_cluster1), mean)) %>%
  pivot_longer(cols = -time_pt, names_to = "gene", values_to = "expression")

# Refine trend classification 
trend_data <- gene_expr_data %>%
  group_by(gene) %>%
  arrange(time_pt) %>%
  summarize(
    trend = if (expression[1] > expression[n()] && all(diff(expression) <= 0)) {
      "Downward Trend"  
    } else if (expression[1] < expression[n()] && all(diff(expression) >= 0)) {
      "Upward Trend"     
    } else if (which.max(expression) %in% 3:(n() - 2) &&
               expression[1] < max(expression) && expression[n()] < max(expression)) {
      "Peaked Trend"    
    } else if (which.min(expression) %in% 3:(n() - 2) &&
               expression[1] > min(expression) && expression[n()] > min(expression)) {
      "Valley Trend"    
    } else if (sd(expression) < 0.2) {
      "Flat Trend"      
    } else {
      "No Clear Trend"   
    }
  )

# Check for empty categories
print(table(trend_data$trend))  

# Merge the trend classification with the expression data
gene_expr_with_trend <- gene_expr_data %>%
  inner_join(trend_data, by = "gene")

# Plot genes with the same trend
plot_genes_by_trend <- function(data, trend_label) {
  ggplot(data %>% filter(trend == trend_label), aes(x = time_pt, y = expression, color = gene, group = gene)) +
    geom_line() +
    geom_point(size = 2) +
    geom_text_repel(aes(label = gene), size = 3, max.overlaps = 500) +  # Increase max.overlaps for dense labeling
    labs(
      title = paste("Gene Expression Over Time -", trend_label),
      x = "Time Point",
      y = NULL,  # Remove y-axis label
      color = "Gene"
    ) +
    theme_minimal() +
    theme(legend.position = "right")
}


# Generate and save plots
plot_upward <- plot_genes_by_trend(gene_expr_with_trend, "Upward Trend")
ggsave(filename = "G:/My Drive/Viv_MACHS/R/upward_trend_genes_cluster1.png", plot = plot_upward, width = 10, height = 8, dpi = 300, bg = "white")

plot_downward <- plot_genes_by_trend(gene_expr_with_trend, "Downward Trend")
ggsave(filename = "G:/My Drive/Viv_MACHS/R/downward_trend_genes_cluster1.png", plot = plot_downward, width = 10, height = 8, dpi = 300, bg = "white")

plot_peaked <- plot_genes_by_trend(gene_expr_with_trend, "Peaked Trend")
ggsave(filename = "G:/My Drive/Viv_MACHS/R/peaked_trend_genes_cluster1.png", plot = plot_peaked, width = 10, height = 8, dpi = 300, bg = "white")

plot_valley <- plot_genes_by_trend(gene_expr_with_trend, "Valley Trend")
ggsave(filename = "G:/My Drive/Viv_MACHS/R/valley_trend_genes_cluster1.png", plot = plot_valley, width = 10, height = 8, dpi = 300, bg = "white")

plot_flat <- plot_genes_by_trend(gene_expr_with_trend, "Flat Trend")
ggsave(filename = "G:/My Drive/Viv_MACHS/R/flat_trend_genes_cluster1.png", plot = plot_flat, width = 10, height = 8, dpi = 300, bg = "white")

plot_no_clear_trend <- plot_genes_by_trend(gene_expr_with_trend, "No Clear Trend")
ggsave(filename = "G:/My Drive/Viv_MACHS/R/no_clear_trend_genes_cluster1.png", plot = plot_no_clear_trend, width = 10, height = 8, dpi = 300, bg = "white")


         ###################### TREND SEPERATION FOR CLUSTER 3 #####################
#Identify the top marker genes for Cluster 3
seurat_integrated_genelist <- PrepSCTFindMarkers(seurat_integrated_genelist)
cluster3_genes <- FindMarkers(seurat_integrated_genelist, ident.1 = 3, logfc.threshold = 0.25)
top_genes_cluster3 <- rownames(cluster3_genes)[1:100]

# Prepare data for line plot
gene_expr_data <- FetchData(seurat_integrated_genelist, vars = c("time_pt", top_genes_cluster3)) %>%
  group_by(time_pt) %>%
  summarize(across(all_of(top_genes_cluster3), mean)) %>%
  pivot_longer(cols = -time_pt, names_to = "gene", values_to = "expression")

# Refine trend classification 
trend_data <- gene_expr_data %>%
  group_by(gene) %>%
  arrange(time_pt) %>%
  summarize(
    trend = if (expression[1] > expression[n()] && all(diff(expression) <= 0)) {
      "Downward Trend"  
    } else if (expression[1] < expression[n()] && all(diff(expression) >= 0)) {
      "Upward Trend"     
    } else if (which.max(expression) %in% 3:(n() - 2) &&
               expression[1] < max(expression) && expression[n()] < max(expression)) {
      "Peaked Trend"    
    } else if (which.min(expression) %in% 3:(n() - 2) &&
               expression[1] > min(expression) && expression[n()] > min(expression)) {
      "Valley Trend"    
    } else if (sd(expression) < 0.2) {
      "Flat Trend"      
    } else {
      "No Clear Trend"   
    }
  )

# Check for empty categories
print(table(trend_data$trend))  

# Merge the trend classification with the expression data
gene_expr_with_trend <- gene_expr_data %>%
  inner_join(trend_data, by = "gene")

# Plot genes with the same trend
plot_genes_by_trend <- function(data, trend_label) {
  ggplot(data %>% filter(trend == trend_label), aes(x = time_pt, y = expression, color = gene, group = gene)) +
    geom_line() +
    geom_point(size = 2) +
    geom_text_repel(aes(label = gene), size = 3, max.overlaps = 500) +  # Increase max.overlaps for dense labeling
    labs(
      title = paste("Gene Expression Over Time -", trend_label),
      x = "Time Point",
      y = NULL,  # Remove y-axis label
      color = "Gene"
    ) +
    theme_minimal() +
    theme(legend.position = "right")
}


# Generate and save plots
plot_upward_cl3 <- plot_genes_by_trend(gene_expr_with_trend, "Upward Trend")
ggsave(filename = "G:/My Drive/Viv_MACHS/R/upward_trend_genes_cluster3.png", plot = plot_upward_cl3, width = 10, height = 8, dpi = 300, bg = "white")

plot_downward_cl3 <- plot_genes_by_trend(gene_expr_with_trend, "Downward Trend")
ggsave(filename = "G:/My Drive/Viv_MACHS/R/downward_trend_genes_cluster3.png", plot = plot_downward_cl3, width = 10, height = 8, dpi = 300, bg = "white")

plot_peaked_cl3 <- plot_genes_by_trend(gene_expr_with_trend, "Peaked Trend")
ggsave(filename = "G:/My Drive/Viv_MACHS/R/peaked_trend_genes_cluster3.png", plot = plot_peaked_cl3, width = 10, height = 8, dpi = 300, bg = "white")

plot_valley_cl3 <- plot_genes_by_trend(gene_expr_with_trend, "Valley Trend")
ggsave(filename = "G:/My Drive/Viv_MACHS/R/valley_trend_genes_cluster3.png", plot = plot_valley_cl3, width = 10, height = 8, dpi = 300, bg = "white")

plot_flat_cl3 <- plot_genes_by_trend(gene_expr_with_trend, "Flat Trend")
ggsave(filename = "G:/My Drive/Viv_MACHS/R/flat_trend_genes_cluster3.png", plot = plot_flat_cl3, width = 10, height = 8, dpi = 300, bg = "white")

plot_no_clear_trend_cl3 <- plot_genes_by_trend(gene_expr_with_trend, "No Clear Trend")
ggsave(filename = "G:/My Drive/Viv_MACHS/R/no_clear_trend_genes_cluster3.png", plot = plot_no_clear_trend_cl3, width = 10, height = 8, dpi = 300, bg = "white")

