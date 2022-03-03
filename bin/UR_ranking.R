# Sandra Lilja
# Built on R version 4.0.4

#' Rank URs
#' 
#' Input directory: 
#' Should contain subdirectories, one for each time point to be included (time point need to be in subdirectory name). 
#' Ensure that the directory  only contains files included in this analysis, or error will occur when listing files.
#' Ensure the names of files follows the system <[cellType]_[timePoint]_[other description].xls>
#' Ensure that the cell types are named equally between the different time points
#' 
#' Output directory:
#' Note that output files which already exists will be overwritten
#' 
#' @param dir.data E.g., data/IPA_UR-prediction
#' @param dir.out E.g., data/UR-rank
#' 
#' @return no return, just saves the files into data folder
#' @export
#'


library(readxl)
library(dplyr)

UR_ranking <- function(dir.data, dir.out){
  
  molecule_types <- c('complex', 'cytokine', 'G-protein coupled receptor',
                      'group', 'growth factor', 'ligand-dependent nuclear receptor',
                      'transmembrane receptor') 
  ipa_columns <- c("Upstream Regulator", "Molecule Type",
                   "Activation z-score", "p-value of overlap")
  
  
  # Rank the URs based on in how many cell types and at how many time points they were predicted in
  
  # read in the ipa output files
  ipa_files <- list.files(dir.data, recursive = T, full.names = T)
  
  ipa <- c()
  for (i in 1:length(ipa_files)){
    ipa[[i]] <- read_excel(ipa_files[i])
    ipa[[i]] <- ipa[[i]][,1:9]
    if (length(grep('All rights reserved', colnames(ipa[[i]])[1]))==1){
      colnames(ipa[[i]]) <- ipa[[i]][1,]
      ipa[[i]] <- ipa[[i]][-1,]
    }
    ipa[[i]] <- ipa[[i]][,which(colnames(ipa[[i]]) %in% ipa_columns)]
    
    # sort the data to only include significant predictions
    if (any(is.na(ipa[[i]]$`Activation z-score`))){
      ipa[[i]]$`Activation z-score`[which(is.na(ipa[[i]]$`Activation z-score`))] <- 0 # NAs are considered not significant
      ipa[[i]]$`Activation z-score` <- as.numeric(ipa[[i]]$`Activation z-score`)
    }
    if (any(is.na(ipa[[i]]$`p-value of overlap`))){
      ipa[[i]]$`p-value of overlap`[which(is.na(ipa[[i]]$`p-value of overlap`))] <- 1 # NAs are considered not significant
      ipa[[i]]$`p-value of overlap` <- as.numeric(ipa[[i]]$`p-value of overlap`)
    }
    
    ipa[[i]] <- ipa[[i]][which(abs(ipa[[i]]$`Activation z-score`) >= 2),]
    ipa[[i]] <- ipa[[i]][which(ipa[[i]]$`p-value of overlap` <= 0.05),]
    # sort the data only to include URs of interest
    ipa[[i]] <- ipa[[i]][which(ipa[[i]]$`Molecule Type` %in% molecule_types),]
  }
  celltypes <- sapply(strsplit(ipa_files, '/'), tail, 1)
  celltypes <- gsub(' ', '', celltypes)
  time_points <- sapply(strsplit(celltypes, '_'), '[[', 2)
  celltypes <- sapply(strsplit(celltypes, '_'), '[[', 1)
  names(ipa) <- paste(celltypes, time_points, sep = '_')
  
  
  # Summarize z-scores for each cellType_timePoint and UR
  for (i in 1:length(ipa)){
    ipa[[i]]$'Sample' <- names(ipa)[i]
  }
  ipa_comb <- rbind(ipa[[1]], ipa[[2]])
  if (length(ipa) > 2){
    for (i in 3:length(ipa)){
      ipa_comb <- rbind(ipa_comb, ipa[[i]])
      
    }
  }
  
  # count the number of occurrences for each UR
  URs <- unique(sort(ipa_comb$`Upstream Regulator`))
  UR_rank <- as.data.frame(matrix(NA, ncol = 2, nrow = length(URs)))
  colnames(UR_rank) <- c('UR', 'N cell types and time points')
  UR_rank$UR <- URs
  for (i in 1:length(UR_rank$UR)){
    UR_rank$`N cell types and time points`[i] <- length(which(ipa_comb$`Upstream Regulator` == UR_rank$UR[i]))
  }
  UR_rank <- UR_rank[order(UR_rank$`N cell types and time points`, decreasing = T),]
  
  print('write ranked list to out')
  # write.csv(ipa_comb, paste(dir.out, '/combined_UR_regulations.csv', sep = ''), row.names = F)
  write.csv(UR_rank, paste(dir.out, '/UR_ranking.csv', sep = ''), row.names = F)
  
}
