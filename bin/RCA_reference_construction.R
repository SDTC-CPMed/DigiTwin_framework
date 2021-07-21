# Sandra Lilja

#######################
### GlobalPanel PBMC data ###
library(preprocessCore)
library(dplyr)
library(AnnotationDbi)
library("hgu133a.db")
rm(list=ls())
# Why 
expr = read.table('/data/sanli71/RCA_reference_playaround/GlobalPanel_RAW_data/U133AGNF1B.gcrma.csv', header = T, row.names = 1, sep = ',')
expr[1:5,1:5]
expr[,length(colnames(expr))+1] <- rownames(expr)
genes <- select(hgu133a.db, expr[,177], c("ENTREZID"))
expr[1:5,175:177]
head(genes)
expr <- full_join(expr, genes, by = c('V177' = 'PROBEID'))
expr[1:5,175:178]
head(sort(expr$ENTREZID))
# test1 <- subset(expr, !duplicated(expr$ENTREZID))
expr <- expr[!(duplicated(expr$ENTREZID) | duplicated(expr$ENTREZID, fromLast =T)),]
rownames(expr) <- expr[,178]
remove(genes)
expr <- expr[,-length(colnames(expr))]
expr <- expr[,-length(colnames(expr))]
# backup <- expr
expr <- as.matrix(expr)
coln <- colnames(expr)
rown <- rownames(expr)
expr <- normalize.quantiles(expr) # Check an extra time so this is not already done on the data
colnames(expr) <- coln
rownames(expr) <- rown
expr[1:5,1:5]
remove(coln, rown)

## Monocytes vs Dendritic vs B vs TNK
XX <- expr
XX[1:5,1:5]
# colnames(XX)[which(grepl('Mono', colnames(XX))==FALSE)] <- 
#   paste('Lympho', colnames(XX)[which(grepl('Mono', colnames(XX))==FALSE)], sep = '_')
colnames(XX)
XX <- XX[,c(grep('Monocytes', colnames(XX)), grep('_B', colnames(XX)), grep('_T', colnames(XX)), grep('NK', colnames(XX)), grep('Dentritic', colnames(XX)))]
colnames(XX)
# XX <- XX[,c(grep('Macrophages', colnames(XX)), grep('_B_lymphoblasts', colnames(XX)), grep('_T', colnames(XX)), grep('NK', colnames(XX)), grep('Dendritic', colnames(XX)))]
colnames(XX)[grep('Tcells', colnames(XX))] <- 
  paste('TNK', colnames(XX)[grep('Tcells', colnames(XX))], sep = '_')
colnames(XX)[grep('NKCells', colnames(XX))] <- 
  paste('TNK', colnames(XX)[grep('NKCells', colnames(XX))], sep = '_')
colnames(XX)[grep('Monocytes', colnames(XX))] <- 
  paste('Monocytes', colnames(XX)[grep('Monocytes', colnames(XX))], sep = '_')
colnames(XX)[grep('_B', colnames(XX))] <- 
  paste('BCells', colnames(XX)[grep('_B', colnames(XX))], sep = '_')
colnames(XX)[grep('Dentritic', colnames(XX))] <- 
  paste('Dendritic', colnames(XX)[grep('Dentritic', colnames(XX))], sep = '_')
colnames(XX)
### Contruct the reference
XX <- sapply(split(seq_len(ncol(XX)), sapply(strsplit(colnames(XX), '_'), '[[', 1)), function(cis) rowMeans(XX[,cis,drop=F]))
mediana = apply(X = XX,MARGIN = 1,FUN = median)
FCs = log10(XX/mediana) # change to log2?
co <- 0.5
colSums(FCs>=co)
IN = rowSums(FCs>=co)!=0 
sum(IN)
Tcellref = log10(XX[IN,])
write.csv(Tcellref, '/data/sanli71/SAR_allergen_challenge_timeseries/SAR_healty_control_study/RCA_out_210211/RCA_references/GlobalPanel_MonvsDCvsBvsTNK.csv')
### Saving a new reference to RCA sysdata
data(sysdata, package='RCA')
sysdata = sysdata[!(names(sysdata) %in% c('GlobalPanel_MonvsDCvsBvsTNK'))]
sysdata = append(sysdata,
                 list('GlobalPanel_MonvsDCvsBvsTNK'=Tcellref))
save(sysdata, file="/home/sanli71/R/x86_64-pc-linux-gnu-library/3.4/RCA/data/sysdata.rda")
names(sysdata)


### Huans PBMC data ###
rm(list=ls())
expr = read.csv('/data/sanli71/RCA_reference_playaround/PBMC_Huan_RAW_data/pbmc_huan_samples/normalized/normalized_matrix_by_entrez-id_with_sample-names.tsv', header = T, row.names = 1, sep = '\t')
expr <- 2^expr # remove log2 transformation from the normalized data
expr <- expr[,-grep('PBMC', colnames(expr))]
# genes = rownames(expr)

## CD8 vs CD4 vs NK cells
XX <- expr
XX <- XX[,-grep('Mono', colnames(XX))]
XX <- XX[,-grep('T', colnames(XX))]
XX <- XX[,-grep('B.cell', colnames(XX))]
colnames(XX)
### Contruct the reference
XX <- sapply(split(seq_len(ncol(XX)), sapply(strsplit(colnames(XX), '_'), '[[', 1)), function(cis) rowMeans(XX[,cis,drop=F]))
mediana = apply(X = XX,MARGIN = 1,FUN = median)
FCs = log10(XX/mediana)
co <- 0.3
colSums(FCs>=co)
IN = rowSums(FCs>=co)!=0 
sum(IN)
Tcellref = log10(XX[IN,])
write.csv(Tcellref, '/data/sanli71/SAR_allergen_challenge_timeseries/SAR_healty_control_study/RCA_out_210211/RCA_references/CellRef_Huan_PBMC_CD4vsCD8vsNK.csv')
### Saving a new reference to RCA sysdata
data(sysdata, package='RCA')
sysdata = sysdata[!(names(sysdata) %in% c('CellRef_Huan_PBMC_CD4vsCD8vsNK'))]
sysdata = append(sysdata,
                 list('CellRef_Huan_PBMC_CD4vsCD8vsNK'=Tcellref))
save(sysdata, file="/home/sanli71/R/x86_64-pc-linux-gnu-library/3.4/RCA/data/sysdata.rda")
names(sysdata)

## NT vs Th1 vs Th2 vs Th17 vs Treg
XX <- expr
XX <- XX[,-grep('Mono', colnames(XX))]
XX <- XX[,-grep('B.cell', colnames(XX))]
XX <- XX[,-grep('CD', colnames(XX))]
XX <- XX[,-grep('NK', colnames(XX))]
colnames(XX)
### Contruct the reference
XX <- sapply(split(seq_len(ncol(XX)), sapply(strsplit(colnames(XX), '_'), '[[', 1)), function(cis) rowMeans(XX[,cis,drop=F]))
mediana = apply(X = XX,MARGIN = 1,FUN = median)
FCs = log10(XX/mediana)
# co <- 0.3
co <- 1
colSums(FCs>=co)
IN = rowSums(FCs>=co)!=0 
sum(IN)
Tcellref = log10(XX[IN,])
write.csv(Tcellref, '/data/sanli71/SAR_allergen_challenge_timeseries/SAR_healty_control_study/RCA_out_210211/RCA_references/CellRef_Huan_PBMC_CD4subsets.csv')
### Saving a new reference to RCA sysdata
data(sysdata, package='RCA')
sysdata = sysdata[!(names(sysdata) %in% c('CellRef_Huan_PBMC_CD4subsets'))]
sysdata = append(sysdata,
                 list('CellRef_Huan_PBMC_CD4subsets'=Tcellref))
save(sysdata, file="/home/sanli71/R/x86_64-pc-linux-gnu-library/3.4/RCA/data/sysdata.rda")
names(sysdata)

