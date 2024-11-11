# Load standard (top 2000 variable genes)
seurat_integrated_standard <- FindNeighbors(seurat_integrated_standard, dims = 1:10, reduction = "harmony")
seurat_integrated_standard <- FindClusters(seurat_integrated_standard, resolution = 0.2, cluster.name = "Standard_clusters")
seurat_integrated_standard <- RunUMAP(seurat_integrated_standard, dims = 1:10, reduction = "harmony")

seurat_integrated_standard$time_pt <- factor(seurat_integrated_standard$time_pt, levels = c("zero", "two", "four", "six", "twenty4", "fourty8"))

st_c<- DimPlot(seurat_integrated_standard, label=T, label.box=T) + NoLegend()
st_condition<- DimPlot(seurat_integrated_standard, group.by = "condition") 
st_time<- DimPlot(seurat_integrated_standard, group.by = "condition_time") 
st_condition_time<- DimPlot(seurat_integrated_standard, group.by = "time_pt") 
st_study<- DimPlot(seurat_integrated_standard, group.by = "study") 

### count###
# counts <- prop.table(table(Idents(seurat_integrated_standard), seurat_integrated_standard$orig.ident), margin = 2)
# #write proportion per cell type
# write.csv(counts, paste0(output_integrated,"_seurat_integrated_standard_counts.csv"), row.names =TRUE)

##########################################################
seurat_integrated_genelist <- FindNeighbors(seurat_integrated_genelist, dims = 1:10, reduction = "harmony")
seurat_integrated_genelist <- FindClusters(seurat_integrated_genelist, resolution = 0.2, cluster.name = "AMDSG_PSG")
seurat_integrated_genelist <- RunUMAP(seurat_integrated_genelist, dims = 1:10, reduction = "harmony")

seurat_integrated_genelist$time_pt <- factor(seurat_integrated_genelist$time_pt, levels = c("zero", "two", "four", "six", "twenty4", "fourty8"))

sg_c<- DimPlot(seurat_integrated_genelist, label=T, label.box=T) + NoLegend()
sg_condition<- DimPlot(seurat_integrated_genelist, group.by = "condition") 
sg_time<- DimPlot(seurat_integrated_genelist, group.by = "condition_time") 
sg_condition_time<- DimPlot(seurat_integrated_genelist, group.by = "time_pt") 
sg_study<- DimPlot(seurat_integrated_genelist, group.by = "study" ) 

### count###
# counts <- prop.table(table(Idents(seurat_integrated_genelist), seurat_integrated_genelist$orig.ident), margin = 2)
# #write proportion per cell type
# write.csv(counts, paste0(output_integrated,"_seurat_integrated_genelist_counts.csv"), row.names =TRUE)
##########################################################

# Count table
condition_counts <- table(Idents(seurat_integrated_genelist), seurat_integrated_genelist$condition)
print(condition_counts)

# Percentage table
condition_percentages <- prop.table(condition_counts, margin = 2) * 100
print(round(condition_percentages, 2))

# Save to CSV
write.csv(condition_counts, "G:/My Drive/Viv_MACHS/R/markers_merged/condition_counts_amdsg_psg.csv")  
write.csv(condition_percentages, "G:/My Drive/Viv_MACHS/R/markers_merged/condition_percentages_amdsg_psg.csv")

## condition_Time
# Count table
condition_time_counts <- table(Idents(seurat_integrated_genelist), seurat_integrated_genelist$condition_time)
print(condition_time_counts)

# Percentage table
condition_time_percentages <- prop.table(condition_time_counts, margin = 2) * 100
print(round(condition_time_percentages, 2))

# Save to CSV
write.csv(condition_time_counts, "G:/My Drive/Viv_MACHS/R/markers_merged/condition_time_counts_amdsg_psg.csv")
write.csv(condition_time_percentages, "G:/My Drive/Viv_MACHS/R/markers_merged/condition_time_percentages_amdsg_psg.csv")

# Count table
time_counts <- table(Idents(seurat_integrated_genelist), seurat_integrated_genelist$time_pt)    
print(time_counts)

# Percentage table
time_percentages <- prop.table(time_counts, margin = 2) * 100
print(round(time_percentages, 2))

# Save to CSV
write.csv(time_counts, "G:/My Drive/Viv_MACHS/R/markers_merged/time_counts_amdsg_psg.csv")
write.csv(time_percentages, "G:/My Drive/Viv_MACHS/R/markers_merged/time_percentages_amdsg_psg.csv")

