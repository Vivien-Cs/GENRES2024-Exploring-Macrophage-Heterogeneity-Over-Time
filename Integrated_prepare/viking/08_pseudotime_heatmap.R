

### load heatmap ###



# Load additional gene lists
cell_extrinsic <- read_excel("G:/My Drive/Genres2024/R/Cell_extrinsic_VC.xlsx")
cell_intrinsic <- read_excel("G:/My Drive/Genres2024/R/Cell_intrinsic_VC.xlsx")
AMDSG_48hr <- read_excel("G:/My Drive/Genres2024/R/AMDSG_48hr.xlsx")
PSG_48hr <- read_excel("G:/My Drive/Genres2024/R/PSG_48hr.xlsx")

# Combine gene lists
combined_genes <- rbind(
  cell_extrinsic %>% mutate(TimeSpecific = "All"),
  cell_intrinsic %>% mutate(TimeSpecific = "All"),
  AMDSG_48hr %>% mutate(TimeSpecific = "48hr"),
  PSG_48hr %>% mutate(TimeSpecific = "48hr")
)

# Filter genes to plot
genes_to_plot <- combined_genes$Gene[combined_genes$Gene %in% rownames(seurat_merged)]
gene_annotation <- combined_genes$Function[combined_genes$Gene %in% rownames(seurat_merged)]
time_specific <- combined_genes$TimeSpecific[combined_genes$Gene %in% rownames(seurat_merged)]

# Extract expression data
expression_data <- as.matrix(GetAssayData(seurat_merged, assay = "RNA", slot = "data")[genes_to_plot, ])

# Ensure pseudotime order is available
pseudotime_order <- order(seurat_merged$time_pt, decreasing = FALSE)

# Create heatmap
heatmap <- Heatmap(
  expression_data[, pseudotime_order],
  name = "Expression",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_title = "Pseudotime",
  row_title = "Genes",
  left_annotation = rowAnnotation(
    Function = gene_annotation,
    TimeSpecific = time_specific,
    col = list(
      Function = setNames(rainbow(length(unique(combined_genes$Function))), unique(combined_genes$Function)),
      TimeSpecific = c("All" = "black", "48hr" = "red")
    ),
    show_legend = TRUE
  ),
  col = colorRamp2(c(min(expression_data), mean(expression_data), max(expression_data)), 
                   c("white", "pink", "red")),
  row_split = gene_annotation,
  row_gap = unit(2, "mm"),
  row_names_gp = gpar(fontsize = 10),
  na_col = "grey",
  heatmap_legend_param = list(
    title = "Expression",
    at = c(0, 2, 4, 6, 8),
    labels = c("0", "2", "4", "6", "8"),
    legend_direction = "horizontal",
    legend_width = unit(4, "cm"),
    title_position = "topcenter"
  )
)

# Draw the heatmap
draw(heatmap, heatmap_legend_side = "bottom")

# Save the heatmap
ggsave("heatmap_results/heatmap.png", plot = grid.grabExpr(draw(heatmap, heatmap_legend_side = "bottom")), device = "png", width = 10, height = 8, units = "in", dpi = 300)

# it worksssss or not (expression at 0 so high???)
draw(heatmap)




                                               ###### filter to run specific functions###############

# Define the selected function to filter genes
selected_function <- "Cytokine_production"

# Load gene lists
cell_extrinsic <- read_excel("G:/My Drive/Genres2024/R/Cell_extrinsic_VC.xlsx")
cell_intrinsic <- read_excel("G:/My Drive/Genres2024/R/Cell_intrinsic_VC.xlsx")
AMDSG_48hr <- read_excel("G:/My Drive/Genres2024/R/AMDSG_48hr.xlsx")
PSG_48hr <- read_excel("G:/My Drive/Genres2024/R/PSG_48hr.xlsx")

# Select relevant columns and combine gene lists
cell_extrinsic <- read_excel("G:/My Drive/Genres2024/R/Cell_extrinsic_VC.xlsx")
cell_intrinsic <- read_excel("G:/My Drive/Genres2024/R/Cell_intrinsic_VC.xlsx")
AMDSG_48hr <- read_excel("G:/My Drive/Genres2024/R/AMDSG_48hr.xlsx")
PSG_48hr <- read_excel("G:/My Drive/Genres2024/R/PSG_48hr.xlsx")

# Combine gene lists
combined_genes <- rbind(
  cell_extrinsic %>% mutate(TimeSpecific = "All"),
  cell_intrinsic %>% mutate(TimeSpecific = "All"),
  AMDSG_48hr %>% mutate(TimeSpecific = "48hr"),
  PSG_48hr %>% mutate(TimeSpecific = "48hr")
)

# Filter genes and annotations based on the selected function
selected_function <- "Cytokine_production"
filtered_genes <- combined_genes %>%
  filter(Function == selected_function) %>%
  select(Gene, Function, TimeSpecific)

# check that filtered genes are present in the Seurat object
genes_to_plot <- filtered_genes$Gene[filtered_genes$Gene %in% rownames(seurat_merged)]
gene_annotation <- filtered_genes$Function[filtered_genes$Gene %in% rownames(seurat_merged)]
time_specific <- filtered_genes$TimeSpecific[filtered_genes$Gene %in% rownames(seurat_merged)]

# Extract expression data for the filtered genes
expression_data <- as.matrix(GetAssayData(seurat_merged, assay = "RNA", slot = "data")[genes_to_plot, ])

# Ensure pseudotime order is available
pseudotime_order <- order(seurat_merged$time_pt, decreasing = FALSE)

## Define color palette for the function and timespecific annotation
function_colors <- setNames(rainbow(length(unique(filtered_genes$Function))), unique(filtered_genes$Function))
time_specific_colors <- c("All" = "black", "48hr" = "blue")

# Create and draw the heatmap
heatmap <- Heatmap(
  expression_data[, pseudotime_order],
  name = "Expression",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_title = "Pseudotime",
  row_title = "Genes",
  left_annotation = rowAnnotation(
    Function = gene_annotation,
    TimeSpecific = time_specific,
    col = list(
      Function = function_colors,
      TimeSpecific = time_specific_colors
    ),
    show_legend = TRUE
  ),
  col = colorRamp2(c(min(expression_data, na.rm = TRUE), 0, mean(expression_data, na.rm = TRUE), max(expression_data, na.rm = TRUE)), 
                   c("white", "white", "pink", "red")),
  row_split = gene_annotation,
  row_gap = unit(2, "mm"),
  row_names_gp = gpar(fontsize = 10),
  na_col = "grey",
  heatmap_legend_param = list(
    title = "Expression",
    at = c(0, 2, 4, 6, 8),
    labels = c("0", "2", "4", "6", "8"),
    legend_direction = "horizontal",
    legend_width = unit(4, "cm"),
    title_position = "topcenter"
  )
)

# Draw the heatmap with adjusted layout
draw(heatmap, heatmap_legend_side = "bottom")

# Save the heatmap as a PNG file in the specified folder yaaay it worked??
ggsave("heatmap_results/heatmap_cytokine_production.png", plot = grid.grabExpr(draw(heatmap, heatmap_legend_side = "bottom")), device = "png", width = 10, height = 8, units = "in", dpi = 300)

                                                 ######################### specific genes ########################

# Define the gene of interest
gene_of_interest <- "Nos2"

# Extract expression data for the gene of interest
expression_data <- FetchData(seurat_merged, vars = gene_of_interest)
print(head(expression_data))

# Combine the expression data with the timepoints and conditions
expression_data <- data.frame(Expression = expression_data[[gene_of_interest]], 
                              TimePoint = seurat_merged@meta.data$time_pt, 
                              Condition = seurat_merged@meta.data$condition)
print(head(expression_data))

# Convert time points to factor with specific order
expression_data$TimePoint <- factor(expression_data$TimePoint, levels = c("zero", "six", "twenty4", "fourty8"))

# Summarise the mean expression for each condition and time point
expression_summary <- expression_data %>%
  group_by(Condition, TimePoint) %>%
  summarise(mean_expression = mean(Expression, na.rm = TRUE),
            sd_expression = sd(Expression, na.rm = TRUE))
print(expression_summary)

# Creating plot
ggplot(expression_summary, aes(x = TimePoint, y = mean_expression, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_expression - sd_expression, ymax = mean_expression + sd_expression), 
                position = position_dodge(width = 0.9), width = 0.25) +
  labs(title = paste("Mean Expression of", gene_of_interest, "Over Time in M0 and M1 Macrophages"), 
       x = "Time Point (hours)", 
       y = "Expression Level",
       fill = "Condition") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#density plt #to be fixed
ggplot(expression_data, aes(x = Expression, fill = Condition, color = Condition)) +
  geom_density(alpha = 0.5, adjust = 1.5) +
  facet_wrap(~ TimePoint, scales = "free_x") +
  labs(title = paste("Density of", gene_of_interest, "Expression Over Time in M0 and M1 Macrophages"), 
       x = "Expression Level", 
       y = "Density",
       fill = "Condition",
       color = "Condition") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
