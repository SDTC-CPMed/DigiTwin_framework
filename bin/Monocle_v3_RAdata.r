# Sandra Lilja

library('monocle')
library('reshape2')
library(R.utils)
rm(list=ls())

# sampleid <- 'HC'
sampleid <- 'HA'

dir.home <- getwd()
dir.data <- paste(dir.home, '/MONOCLE/data/DEGS_with_Monocle_', sampleid, '/Matrix_in/', sep = '')
dir.out <- paste(dir.home, '/MONOCLE/data/DEGS_with_Monocle_', sampleid, '/Monocle_out/', sep = '')

## load data:
scrna = read.csv(paste(dir.data, sampleid, '_ENTREZ_expression_matrix.csv', sep = ''), sep = ' ', row.names = 1)
scrna[1:5,1:5]

# Note:
# Dependent on if the calculations are between healthy and sick, or allergen challenged and diluent, or anything else, 
# remove the cells which are not included

### SAR remove timepoint 0 - for allergen challenged vs diluent challenged calculations
scrna <- scrna[,!grepl('_0h_C_', colnames(scrna))]

## Prepare groups information
pd = data.frame(keys=colnames(scrna))
pd$keys2 = sapply(strsplit(colnames(scrna),'__'),function(Z) Z[1])
pd$kosz = sapply(strsplit(pd$keys2,'_'),function(Z) Z[1])
pd$sampleSource =  sapply(strsplit(pd$keys2,'_'), '[', 3) # eg tissue or timepoint
pd$celltype = sapply(pd$kosz, function(x) substr(x,1,nchar(x))) # eg CD4
pd$state = sapply(strsplit(pd$keys2,'_'), '[', 4) # eg RA vs healthy
pd$sampleId = sapply(strsplit(pd$keys2,'_'),function(Z) Z[length(Z)]) # cellular barcode
pd$kosz = NULL
rownames(pd) = colnames(scrna)

### ### #### ### ###
## START MONOCLE: ##
### ### #### ### ###

## create variables:
fd = data.frame(gene_short_name = rownames(scrna))
rownames(fd) = rownames(scrna)
pd <- new("AnnotatedDataFrame", data = pd)
fd <- new("AnnotatedDataFrame", data = fd)
# Recommended for UMI count data
scrna <- newCellDataSet(as(as.matrix(scrna),'sparseMatrix'),
                        phenoData = pd,
                        featureData = fd,
                        lowerDetectionLimit = 0.5,
                        expressionFamily=negbinomial()) # Negbinomial() is slow so for larger data sets negbinomial.size() is recommended

f_estimateDispersions <- function(x){
  file <- try(estimateDispersions(x))
  if (class(file) == "try-error"){
    print(paste(ucelltypes, usource, 'error during estimateDispersion, try withot this step', sep = ' '))
    file <- x
  }
  return(file)
}

## GET DEGs:
## unchallenged vs challenged
ucelltypes = unique(pData(scrna)$celltype)
ucelltypes = ucelltypes[c(1:4,6:10,5)] ### NOTE: do Th2, Th17 and NK last thing
usource = unique(pData(scrna)$sampleSource)
for (i in 1:length(usource)){
  for (z  in 1 : length(ucelltypes)){
    if (length(grep('TRUE', pData(scrna)$celltype == ucelltypes[z] & 
                    pData(scrna)$sampleSource == usource[i] & 
                    pData(scrna)$state == 'D')) <= 2 | 
        length(grep('TRUE', pData(scrna)$celltype == ucelltypes[z] & 
                    pData(scrna)$sampleSource == usource[i] & 
                    pData(scrna)$state == 'A')) <= 2){
      print(paste(ucelltypes[z], usource[i], "sample size <= 2 in one or both states", sep = ' '))
      next
    }
    scrna_subset = scrna[,pData(scrna)$celltype == ucelltypes[z]]
    scrna_subset = scrna_subset[,scrna_subset@phenoData@data$sampleSource == usource[i]]
    scrna_subset <- estimateSizeFactors(scrna_subset)
    scrna_subset <- f_estimateDispersions(scrna_subset)
    scrna_subset <- detectGenes(scrna_subset, min_expr = 0.1)
    expressed_genes <- row.names(subset(fData(scrna_subset), num_cells_expressed >= 3))
    scrna_subset = scrna_subset[expressed_genes,]

    diff_test_res <- differentialGeneTest(scrna_subset,
                                          fullModelFormulaStr = "~state", cores = 30)

    write.table(diff_test_res,paste(dir.out,'Monocle_DEGs_', sampleid, '_',ucelltypes[z],'_',usource[i],'_AllergenChallenged_vs_NonChallenged.txt',sep=''),row.names=F)
    rm(list=c('scrna_subset','diff_test_res'))
  }
}

## gzip out
filename <- list.files(dir.out, pattern = '.txt', full.names = T)
filename <- filename[!grepl('.gz', filename)]
for (i in 1:length(filename)){
  gzip(filename[i])
}

