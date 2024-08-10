## volcano plots ###
# Find markers between M0 and M1 groups (just for test)
de_genes <- FindMarkers(seurat_merged, 
                        ident.1 = M0_group, 
                        ident.2 = M1_group,
                        only.pos = FALSE,   # Include both up and down-regulated genes if only.pos=FALSE
                        min.pct = 0.25, 
                        logfc.threshold = 0.25)
# p=>0.05
de_genes$p_val_adj[de_genes$p_val_adj == 0] <- 0.05

# Create a volcano plot with adjusted parameters
EnhancedVolcano(de_genes,
                lab = rownames(de_genes),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                selectLab = top20,
                pCutoff = 0.05,
                FCcutoff = 0.25,
                pointSize = 1.0,
                labSize = 3.0,
                max.overlaps = 30,
                drawConnectors = TRUE,
                widthConnectors = 0.25,
                colConnectors = 'grey50',
                title = 'M0 vs M1 Differential Expression',
                subtitle = paste0('Top 20 genes labeled'),
                caption = paste0('Total = ', nrow(de_genes), ' genes'),
                xlim = c(-max(abs(de_genes$avg_log2FC)), max(abs(de_genes$avg_log2FC))),
                ylim = c(0, max(-log10(de_genes$p_val_adj))),
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 4.0)