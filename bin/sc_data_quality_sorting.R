# Sandra Lilja
# R version 3.4

library(Seurat)

sc_data_quality_sorting <- function(filename){
  
  #filename = '/data/sharedData/SAR_allergen_challenge_timeseries/SAR_healty_control_study/DGE_data/HC_min200genesPerCell_UMI_expression_matrix.txt.gz'
  #filenames = c('/data/sharedData/SAR_allergen_challenge_timeseries/SAR_allergic_patient_study/DGE_data/HA_min200genesPerCell_UMI_expression_matrix.txt.gz',
  #              '/data/sharedData/SAR_allergen_challenge_timeseries/SAR_healty_control_study/DGE_data/HC_min200genesPerCell_UMI_expression_matrix.txt.gz')
  
  X <- read.table(filename, header = T, sep = '\t', row.names = 1)
  X[1:5,1:5]
  
  #X.seu <- CreateSeuratObject(X2, min.cells = 3)
  ### Seurat QC and data sorting
  #X.seu <- CreateSeuratObject(X, min.cells = 3, min.genes = 200)
  X.seu <- CreateSeuratObject(X, min.cells = 3, min.features = 200)
  
  
  X.seu[["percent.mt"]] <- PercentageFeatureSet(X.seu, pattern = "^MT-")
#  VlnPlot(X.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
#  plot1 <- FeatureScatter(X.seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
#  plot2 <- FeatureScatter(X.seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#  plot1 + plot2
  
  
  ### Sort the data based on results from previous graph
  #X.seu <- subset(X.seu, nFeature_RNA < 15000)
  
  #X.seu <- subset(X.seu, subset = nFeature_RNA < 16000 & percent.mt < 20 & nCount_RNA > 200 & nFeature_RNA > 400)
  X.seu <- subset(X.seu, subset = nCount_RNA < 16000)
  
  
  
  ### Extract expression matrix
  X <- as.data.frame(X.seu@assays$RNA@counts)
  X[1:5,1:5]
  
  # write.table(X, paste(tissue, '_min200genesPerCell_sorted_expression_matrix.txt', sep = ''),
  #             row.names = T, col.names = T, sep = '\t', quote = F)
  # 
  # write.csv(X, paste(tissue, '_min200genesPerCell_sorted_expression_matrix.csv', sep = ''), quote = F)

  
  return(X)
  
}



# for(tissue in c("HC")){
#   
#   if(tissue == "HA"){
#     filename = '/data/sharedData/SAR_allergen_challenge_timeseries/SAR_allergic_patient_study/DGE_data/HA_min200genesPerCell_UMI_expression_matrix.txt.gz' 
#   }
#   else{
#     filename = '/data/sharedData/SAR_allergen_challenge_timeseries/SAR_healthy_control_study/DGE_data/HC_min200genesPerCell_UMI_expression_matrix.txt.gz' 
#   }
#   
#   X <- read.table(filename, header = T, sep = '\t', row.names = 1)
#   X[1:5,1:5]
#   
#   ### Seurat QC and data sorting
#   X.seu <- CreateSeuratObject(X, min.cells = 3, min.genes = 100)
#   
#   X.seu[["percent.mt"]] <- PercentageFeatureSet(X.seu, pattern = "^MT-")
#   VlnPlot(X.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#   
#   plot1 <- FeatureScatter(X.seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
#   plot2 <- FeatureScatter(X.seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#   plot1 + plot2
#   
#   
#   ### Sort the data based on results from previous graph
#   print(tissue)
#   X.seu <- subset(X.seu, nFeature_RNA < 16000)
#   
#   
#   ### Extract expression matrix
#   X <- as.data.frame(X.seu@assays$RNA@counts)
#   X[1:5,1:5]
#   
#   write.table(X, paste(tissue, '_min200genesPerCell_sorted_expression_matrix.txt', sep = ''),
#               row.names = T, col.names = T, sep = '\t', quote = F)
#   
#   write.csv(X, paste(tissue, '_min200genesPerCell_sorted_expression_matrix.csv', sep = ''), quote = F)
# }





