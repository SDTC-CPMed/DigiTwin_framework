# Sandra Lilja

library(Seurat)
library(dplyr)
library(Matrix)
library(R.utils)

rm(list=ls())
dir.home <- getwd()
dir.home <- paste(dir.home, '/DGE_data/', sep = '')
dir.data <- dir.home
dir.out <- paste(dir.home, 'sorted_DGEs/', sep = '')

tissue <- "HA" # "HC"
projectname <- sapply(strsplit(dir.home, '/'), '[[', 4)
list.files(dir.data, pattern = tissue)

X <- read.table(paste(dir.data, list.files(dir.data, pattern = tissue), sep = ''), header = T, sep = '\t', row.names = 1)
X[1:5,1:5]

### Seurat QC and data sorting
X.seu <- CreateSeuratObject(X, min.cells = 3, min.genes = 100)
VlnPlot(X.seu, c('nGene', 'nUMI'), nCol = 2, size.x.use = 8, size.y.use = 8, size.title.use = 8)
mito.genes <- grep('^MT-', rownames(X), value = T)
percent.mito <- Matrix::colSums(X.seu@raw.data[mito.genes, ], na.rm = T)/Matrix::colSums(X.seu@raw.data, na.rm = T)
X.seu <- AddMetaData(X.seu, percent.mito, 'percent.mito')
VlnPlot(X.seu, c('nGene', 'nUMI', 'percent.mito'), nCol = 3, size.x.use = 8, size.y.use = 8, size.title.use = 8)
par(mfrow = c(1, 2))
GenePlot(X.seu, 'nUMI', 'percent.mito')
GenePlot(X.seu, 'nUMI', 'nGene')
dev.off()
vv <- X.seu@meta.data$nGene/X.seu@meta.data$nUMI
plot(density(vv))

### Sort the data based on results from previous graph
if (tissue == 'HA'){
  print(tissue)
  X.seu <- FilterCells(X.seu, subset.names = 'nUMI',
                       high.thresholds = 15000)
}else {
  print('tissue not defined')
}
X.seu@data[1:5,1:5]

### Extract expression matrix
X <- as.data.frame(as.matrix(X.seu@data))
X[1:5,1:5]

write.table(X, paste(dir.out, paste(tissue, '_min200genesPerCell_sorted_expression_matrix.txt', sep = ''), sep = '/'),
            row.names = T, col.names = T, sep = '\t', quote = F)
filename <- list.files(dir.out, pattern = '.txt', full.names = T)
gzip(filename)
write.csv(X, paste(dir.out,  paste(tissue, '_min200genesPerCell_sorted_expression_matrix.csv', sep = ''), sep = '/'), quote = F)

### QC post-sorting
if (length(grep('qc_results.csv', list.files(dir.home))) == 1){
  QC <- as.matrix(read.csv(paste(dir.home, 'qc_results.csv', sep = '/'), sep = ',', header = T))
  QC <- rbind(QC, NA)
} else {
  QC <- matrix(data = NA, ncol = 11, nrow = 1)
  colnames(QC) <- c('project.name', 'date', 'tissue', 'Ncells', 'Ncells above 500 genes', 'Ncells above 1000 genes', 'mean Ngenes per cell', 'min Ngenes per cell', 'max Ngenes per cell', 'max NUMI per cell', 'mean NUMI per cell')
}

head(X.seu@meta.data)
QC[nrow(QC),1] <- projectname
QC[nrow(QC),2] <- as.character(Sys.Date())
QC[nrow(QC),3] <- tissue
QC[nrow(QC),4] <- length(X.seu@meta.data[,1])
QC[nrow(QC),5] <- length(which(X.seu@meta.data[,1]>500)==TRUE)
QC[nrow(QC),6] <- length(which(X.seu@meta.data[,1]>1000)==TRUE)
QC[nrow(QC),7] <- mean(X.seu@meta.data[,1]) # meanGene
QC[nrow(QC),8] <- min(X.seu@meta.data[,1]) # minGene
QC[nrow(QC),9] <- max(X.seu@meta.data[,1]) # maxGene
QC[nrow(QC),10] <- max(X.seu@meta.data[,2]) # maxUMI
QC[nrow(QC),11] <- mean(X.seu@meta.data[,2]) # meanUMI

write.csv(QC, paste(dir.home, 'qc_results.csv', sep = '/'), row.names = F)

