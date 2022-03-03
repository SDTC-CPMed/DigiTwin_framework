# Sandra Lilja

#' Calculate DEGs
#' 
#' @param tissue Either 'HA' for allergic patients of 'HC' for healthy controls
#' 
#' @return no return, just saves the files into data folder
#' @export
#'


library('monocle')
library('reshape2')
library(R.utils)
library(dplyr)


Monocle_v3_RAdata = function(type){
  dir.data <- 'data/DEGS_with_Monocle/Matrix_in/'
  dir.out <- 'data/DEGS_with_Monocle/Monocle_out/'
  
  ## load data:
  if(type == 'HA_vs_HC'){
    infiles <- list.files(dir.data)
    print('aaa')
    HA_in <- read.csv(paste(dir.data, 'HA_ENTREZ_expression_matrix.csv', sep = ''), sep = ' ', row.names = 1)
    HC_in <- read.csv(paste(dir.data, 'HC_ENTREZ_expression_matrix.csv', sep = ''), sep = ' ', row.names = 1)
    print('bbb')
    colnames(HA_in) <- paste('HA', colnames(HA_in), sep = '_')
    colnames(HC_in) <- paste('HC', colnames(HC_in), sep = '_')
    HA_in$'celltypes' <- rownames(HA_in)
    HC_in$'celltypes' <- rownames(HC_in)
    scrna = full_join(HA_in, HC_in)
    remove(HA_in, HC_in, infiles)
    rownames(scrna) <- scrna$celltypes
    scrna <- scrna[,-grep('celltypes', colnames(scrna))]
    scrna[1:5,1:5]
    
    ## SAR remove NonChallenged cells
    scrna <- scrna[,!grepl('_D_', colnames(scrna))]
    unique(sort(sapply(strsplit(colnames(scrna), '_'), '[[', 4)))
    colnames(scrna) <- gsub('_day1_', '_1D_', colnames(scrna))
    colnames(scrna) <- gsub('_day2_', '_2D_', colnames(scrna))
    colnames(scrna) <- gsub('_day3_', '_3D_', colnames(scrna))
    colnames(scrna) <- gsub('_day5_', '_5D_', colnames(scrna))
    colnames(scrna) <- gsub('_day7_', '_7D_', colnames(scrna))
    unique(sort(sapply(strsplit(colnames(scrna), '_'), '[[', 2)))
    
    ## write matrix to out
    write.table(scrna,paste(dir.data, type, '_ENTREZ_expression_matrix.csv' ,sep=''), sep = ' ', col.names = NA, quote = FALSE)
  } else{
    scrna = read.csv(paste(dir.data, type, '_ENTREZ_expression_matrix.csv', sep = ''), sep = ' ', row.names = 1)
  }
  
  print('loaded')
  

  scrna[1:5,1:5]
  
  # Note:
  # Dependent on if the calculations are between healthy and sick, or allergen challenged and diluent, or anything else, 
  # remove the cells which are not included
  
  ### SAR remove timepoint 0 - for allergen challenged vs diluent challenged calculations
  if(type != 'HA_vs_HC'){
    scrna <- scrna[,!grepl('_0h_C_', colnames(scrna))]
  }
  
  if(type == 'HA_vs_HC'){
     ## for SAR HC vs HA data
     pd = data.frame(keys=colnames(scrna))
     pd$keys2 = sapply(strsplit(colnames(scrna),'__'),function(Z) Z[1])
     pd$sampleSource =  sapply(strsplit(pd$keys2,'_'), '[', 4) # eg tissue or timepoint
     pd$celltype = sapply(strsplit(pd$keys2,'_'),function(Z) Z[2])
     sapply(pd$kosz, function(x) substr(x,1,nchar(x))) # eg CD4
     pd$state = sapply(strsplit(pd$keys2,'_'), '[', 1) # eg RA vs healthy
     rownames(pd) = colnames(scrna)
  } else {
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
  }
  
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
    for (z  in 1:length(ucelltypes)){
      if (length(grep('TRUE', pData(scrna)$celltype == ucelltypes[z] &  
                      pData(scrna)$sampleSource == usource[i] & 
                      pData(scrna)$state == unique(pData(scrna)$state)[1])) <= 2 |  
          length(grep('TRUE', pData(scrna)$celltype == ucelltypes[z] &  
                      pData(scrna)$sampleSource == usource[i] & 
                      pData(scrna)$state == unique(pData(scrna)$state)[2])) <= 2){ 
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
  
      write.table(diff_test_res,paste(dir.out,'Monocle_DEGs_', type, '_',ucelltypes[z],'_',usource[i],'_AllergenChallenged_vs_NonChallenged.txt',sep=''),row.names=F)

      rm(list=c('scrna_subset','diff_test_res')) 
    }
  }
  
  
  ## gzip out
  filename <- list.files(dir.out, pattern = '.txt', full.names = T) 
  filename <- filename[!grepl('.gz', filename)] 
  for (i in 1:length(filename)){
   gzip(filename[i])
  }
}
