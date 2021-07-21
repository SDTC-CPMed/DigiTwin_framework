# Sandra Lilja

### Create CellCount per timepoint/sample table
### Create csv files for FC analysis

rm(list=ls())

# sampleid <- 'HC'
# kval = 21
sampleid <- 'HA'
kval = 14

dir.home <- getwd()
dir.data <- paste(dir.home, '/DGE_data/RCA_out/matrices_withCellTypes/full_matrix/', sep = '')
dir.out1 <- paste(dir.home, '/MONOCLE/data/DEGS_with_Monocle_', sampleid, '/', sep = '')
dir.out2 <- paste(dir.out1, 'Matrix_in/', sep = '')

## load data:
X = read.table(paste(dir.data, sampleid, '_min200genesPerCell_sorted_ENTREZ_expression_matrix.knn-smooth_k', kval, '_und-rm_withCellTypes.txt.gz', sep = ''),
               header = T, row.names = 1)
X[1:5,1:5]

#### CellCounts_per_Sample matrices
#########################################
#### detailed cell count
aa <- unique(sapply(strsplit(colnames(X), '_'), '[[', 1))
bb <- unique(paste(sapply(strsplit(colnames(X), '_'), '[[', 2),
                   sapply(strsplit(colnames(X), '_'), '[[', 3),
                   sapply(strsplit(colnames(X), '_'), '[[', 4), sep = '_'))
MM <- matrix(data = NA, nrow = length(aa), ncol = length(bb))
rownames(MM) <- aa
colnames(MM) <- bb

for (i in 1:length(aa)){
  for (z in 1:length(bb))
    MM[i,z] <- length(grep(aa[i], colnames(X[grep(bb[z], colnames(X))])))
}
MM[1:5,1:5]

write.table(MM, paste(dir.out1, 'CellCounts_per_sample_', sampleid, '_k', kval, '.xls', sep = ''), col.names = NA, row.names = T, quote = F, sep = '\t')

#### summarized cell count
bb <- unique(paste(sapply(strsplit(colnames(X), '_'), '[[', 3),
                   sapply(strsplit(colnames(X), '_'), '[[', 4), sep = '_'))

MM <- matrix(data = NA, nrow = length(aa), ncol = length(bb))
rownames(MM) <- aa
colnames(MM) <- bb

for (i in 1:length(aa)){
  for (z in 1:length(bb))
    MM[i,z] <- length(grep(aa[i], colnames(X[grep(bb[z], colnames(X))])))
}
MM[1:5,1:5]
write.table(MM, paste(dir.out1, 'CellCounts_summary_per_sample_', sampleid, '_k', kval, '.xls', sep = ''), col.names = NA, row.names = T, quote = F, sep = '\t')
#########################################

#### DEG/FC csv file
#########################################
write.table(X, paste(dir.out2, sampleid, '_ENTREZ_expression_matrix.csv', sep = ''), sep = ' ', row.names = T, col.names = NA, quote = F)
#########################################

