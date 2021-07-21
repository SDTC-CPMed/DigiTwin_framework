# Sandra Lilja

library(Seurat)
library(dplyr)
library(Matrix)
library(plyr)
library(R.utils)

rm(list=ls())
dir.home <- getwd()
dir.data <- paste(dir.home, '/scData/', sep = '')
dir.out <- paste(dir.home, '/DGE_data/', sep = '')

tissue2 <- 'HA' # 'HC'


dirs <- list.files(path = dir.data, full.names = T, pattern = 'NS500340')
dirs <- list.files(path = dirs, full.names = T, pattern = '_SAR_')


X <- c(1:length(dirs))
for (i in 1:length(dirs)){
  X[i] <- lapply(paste(dirs[i], '/umi.dge.txt.gz', sep = ''), read.table, header = T, row.names = 1, sep = '\t')
}
#X[[1]][1:5,1:5]
#X[[6]][1:5,1:5]
#X[[39]][1:5,1:5]

sampnames <- sapply(strsplit(dirs, '/'), tail, 1)
sampnames <- paste(sapply(strsplit(sampnames, '_'), '[[', 1),
                   sapply(strsplit(sampnames, '_'), '[[', 2),
                   sapply(strsplit(sampnames, '_'), '[[', 3), sep = '_')
length(sampnames) == length(unique(sort(sampnames)))
length(sampnames) == length(X)


for (i in 1:length(sampnames)){
  colnames(X[[i]]) <- paste(sampnames[i], colnames(X[[i]]), sep = '_')
}
#X[[1]][1:5,1:5]
#X[[5]][1:5,1:5]
#X[[9]][1:5,1:5]

as.data.frame(lapply(X, dim))

# Subset the data based on defined quality cut-offs
for (i in 1:length(X)){
  X.seu <- CreateSeuratObject(X[[i]], min.cells = 3, min.genes = 3)
  mito.genes <- grep('^MT-', rownames(X[[i]]), value = T)
  percent.mito <- Matrix::colSums(X.seu@raw.data[mito.genes, ], na.rm = T)/Matrix::colSums(X.seu@raw.data, na.rm = T)
  X.seu <- AddMetaData(X.seu, percent.mito, 'percent.mito')
  X.seu <- FilterCells(X.seu, subset.names = c('nGene', 'percent.mito', 'nUMI'),
                       low.thresholds = c(200, -Inf, 400), high.thresholds = c(Inf, 0.2, Inf))
  X[[i]] <- X[[i]][,c(colnames(X.seu@data))]
}

as.data.frame(lapply(X, dim))  

head(colnames(X[[1]]))
head(rownames(X[[1]]))
ncells <- sum(as.data.frame(lapply(X, dim))[2,])  

# combine tha matrices from all samples
for (i in 1:length(X)){
  X[[i]]$rn <- rownames(X[[i]])
}
X[[5]][1:5,(length(X[[5]])-5):length(X[[5]])]

X <- join_all(X, by = 'rn', type = 'full')
rownames(X) <- X$rn
X[1:5,1:5]
X$rn <- NULL

length(colnames(X)) == ncells

if (any(is.na(X)) == TRUE){
  X[is.na(X)] <- 0
} else {
  print(any(is.na(X)))
}

# output combined matrix
write.table(X, paste(dir.out, tissue2, '_min200genesPerCell_UMI_expression_matrix.txt', sep = ''), col.names = T, row.names = T, sep = '\t', quote = F)
filename <- list.files(dir.out, pattern = tissue2, full.names = T)
gzip(filename)
