# Sandra Lilja


#' Create CellCount per timepoint/sample table
#' Create csv files for FC analysis
#' 
#' @param tissue Either 'HA' for allergic patients of 'HC' for healthy controls
#' 
#' @return no return, just saves the files into data folder
#' @export
#'


pre_DEG_analysis = function(tissue){
  if(tissue == 'HC'){kval = '21'
  } else{kval = '14'}
  

  dir.data <- 'data/RCA_out/full_matrix/'
  dir.out1 <- 'data/DEGS_with_Monocle/'
  dir.out2 <- 'data/DEGS_with_Monocle/Matrix_in/'
  
  
  
  ## load data:
  X = readRDS(paste(dir.data, tissue, '_min200genesPerCell_sorted_ENTREZ_expression_matrix.knn-smooth_k', kval, '_und-rm_withCellTypes.rds', sep = ''))
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
  
  write.table(MM, paste(dir.out1, 'CellCounts_per_sample_', tissue, '_k', kval, '.xls', sep = ''), col.names = NA, row.names = T, quote = F, sep = '\t')
  
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
  write.table(MM, paste(dir.out1, 'CellCounts_summary_per_sample_', tissue, '_k', kval, '.xls', sep = ''), col.names = NA, row.names = T, quote = F, sep = '\t')
  #########################################
  
  #### DEG/FC csv file  
  #########################################
  write.table(X, paste(dir.out2, tissue, '_ENTREZ_expression_matrix.csv', sep = ''), sep = ' ', row.names = T, col.names = NA, quote = F)

  
  #########################################
}
