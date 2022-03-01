# Sandra Lilja
# R version 3.4

#' Check the quality of the data and remove outliers
#' 
#' @param filename The filename of the txt.gz input file with cells in columns and genes in rows
#' 
#' @return The expression matrix without outliers
#' @export
#'

library(Seurat)

sc_data_quality_sorting <- function(filename){
  
  X <- read.table(filename, header = T, sep = '\t', row.names = 1)
  X[1:5,1:5]
  
  X.seu <- CreateSeuratObject(X, min.cells = 3, min.features = 200)
  
  
  X.seu[["percent.mt"]] <- PercentageFeatureSet(X.seu, pattern = "^MT-")
  
  ### If plots are needed to consider the parameters
  #  VlnPlot(X.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  #  plot1 <- FeatureScatter(X.seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
  #  plot2 <- FeatureScatter(X.seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  #  plot1 + plot2
  
  
  ### Sort the data based on results from previous graph
  #X.seu <- subset(X.seu, nFeature_RNA < 15000)
  
  X.seu <- subset(X.seu, subset = nCount_RNA < 16000 & percent.mt < 20)

  ### Extract expression matrix
  X <- as.data.frame(X.seu@assays$RNA@counts)
  X[1:5,1:5]
  
  return(X)
}







