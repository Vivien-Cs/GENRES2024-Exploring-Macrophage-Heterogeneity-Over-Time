
#Initialise an empty list to store Seurat objects for each sample.
seuratObjList<-list()
#######DO NOT RUN IF THE FOLLOWING IS NOT IN PLACE##############
# Please make sure each of the individual sample RAW matrix, genes and barcodes files
# are available under exp_path 

# Define the samples, corresponding time points, and conditions.
lstSamples <- c("M0_0h","M0_6h","M0_24h","M1_6h","M1_24h")
lstTime <- c("zero","six","twenty4","six","twenty4")
lstConditions <- c("M0","M0","M0","M1","M1")

#Loop through each sample to load the data, create Seurat objects, and perform QC and variable feature selection.
for (sample in 1:length(lstSamples)){
  seuratObj <- Read10X(data.dir = paste0(exp_path, lstConditions[[sample]] ,"/", lstTime[[sample]], "/" ))
  seuratObj <- CreateSeuratObject(counts = seuratObj, project = lstSamples[[sample]], min.cells = 3, min.features = 200)
  seuratObj$condition <- lstConditions[[sample]]
  seuratObj$time_pt <- lstTime[[sample]]
  
# Identify the 2,000 most variable features for downstream analysis.
  seuratObj <- FindVariableFeatures(seuratObj, 
                                    selection.method = "vst", 
                                    nfeatures = 2000)
  
# Visualise feature distributions with violin plots. 
  p1 <- VlnPlot(seuratObj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0.1)
  print(seuratObj)
  p2 <- FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  print(p2+p1)
  #seuratObj <- subset(seuratObj, subset = nFeature_RNA > 500 & nFeature_RNA < 4500 & nCount_RNA < 20000 & nCount_RNA > 1000)
  
  seuratObjList[[lstSamples[[sample]]]] <- seuratObj
  print(seuratObjList[[lstSamples[[sample]]]])
}

#######    load 48hr  ################
# Load 48-hour datasets correctly including M0 at 0 hour
lstSamples_48h <- c("M0_0h_48", "M1_48h")  
lstTime_48h <- c("zero", "fourty8")
lstConditions_48h <- c("M0", "M1")
path_48hr <- paste0(base_path,"/data_raw/GSE117176_48hr_RAW")

for (sample in 1:length(lstSamples_48h)) {
  condition <- lstConditions_48h[[sample]]
  data_dir <- file.path(path_48hr, condition, lstTime_48h[[sample]])
  
  seuratObj <- Read10X(data.dir = data_dir)
  seuratObj <- CreateSeuratObject(counts = seuratObj, project = lstSamples_48h[[sample]], min.cells = 3, min.features = 200)
  seuratObj$condition <- condition
  seuratObj$time_pt <- lstTime_48h[[sample]]
  seuratObj <- FindVariableFeatures(seuratObj, selection.method = "vst", nfeatures = 2000)
  
  p1 <- VlnPlot(seuratObj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0.1)
  print(seuratObj)
  p2 <- FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(p2 + p1)
  
  seuratObjList[[lstSamples_48h[[sample]]]] <- seuratObj
  print(seuratObjList[[lstSamples_48h[[sample]]]])
}

#### load 2,4,6hr dataset ######
lstSamples_2_4_6h <- c("M1_2h_rep1", "M1_4h_rep1", "M1_6h_rep1",
                       "M1_2h_rep2", "M1_4h_rep2", "M1_6h_rep2",
                       "M1_2h_rep3", "M1_4h_rep3", "M1_6h_rep3",
                       "M0_0h_rep1", "M0_0h_rep2", "M0_0h_rep3")
lstConditions_2_4_6h <- c("M1", "M1", "M1", "M1", "M1", "M1", "M1", "M1", "M1","M0","M0","M0")
path_2_4_6hr <- paste0(base_path,"/data_raw/E_MTAB_6754_2_4_6hr_RAW/")
lstTime_2_4_6hr <- c("two", "four", "six","two", "four", "six","two", "four", "six","zero", "zero", "zero")

print(paste("Loading data from: ", data_dir))

for (sample in 1:length(lstSamples_2_4_6h)) {
  condition <- lstConditions_2_4_6h[[sample]]
  data_dir <- file.path(path_2_4_6hr,lstConditions_2_4_6h[[sample]],lstSamples_2_4_6h[[sample]])
  seuratObj <- Read10X(data.dir = data_dir)
  seuratObj <- CreateSeuratObject(counts = seuratObj, project = lstSamples_2_4_6h[[sample]], min.cells = 3, min.features = 200)
  seuratObj$condition <- condition
  seuratObj$time_pt <-lstTime_2_4_6hr[[sample]]
  seuratObj <- FindVariableFeatures(seuratObj, selection.method = "vst", nfeatures = 2000)
  
  p1 <- VlnPlot(seuratObj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0.1)
  print(seuratObj)
  p2 <- FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(p2 + p1)
  
  seuratObjList[[lstSamples_2_4_6h[[sample]]]] <- seuratObj
  print(seuratObjList[[lstSamples_2_4_6h[[sample]]]])
}
