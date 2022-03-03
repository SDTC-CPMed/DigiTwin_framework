# Sandra Lilja

#' FC calculations
#' 
#' @param tissue Either 'HA' for allergic patients of 'HC' for healthy controls
#' 
#' @return no return, just saves the files into data folder
#' @export
#'

library(pscl)
library(foreach)
library(doParallel)
library(dplyr)


sc_FC_zero_infl_neg_binomial = function(type){
  mat.data <- 'data/DEGS_with_Monocle/Matrix_in/'
  dir.data <- 'data/DEGS_with_Monocle/Monocle_out/'
  dir.out <- 'data/DEGS_with_Monocle/Monocle_out_withFCs/'
  
  ### load the data
  X <- read.csv(paste(mat.data, type, '_ENTREZ_expression_matrix.csv', sep = ''), sep = ' ', row.names = 1)
  X[1:5,1:5]
  
  
  ### SAR remove timepoint 0
  if(type != 'HA_vs_HC'){
    X <- X[,!grepl('_0h_C_', colnames(X))]
  }
  # ### SAR remove unchallenged group D
  # X <- X[,!grepl('_D_', colnames(X))]
  
  ### cell types
  clls <- unique(sapply(strsplit(colnames(X), '_'), '[[', 1))
  timepoints <- unique(sapply(strsplit(colnames(X), '_'), '[[', 3))
  
  ######### define functions #########
  f_deg <- function(cllx, timepointx){
    nfit <- paste('Monocle_DEGs_', type, '_', cllx, '_', timepointx, 
                  '_AllergenChallenged_vs_NonChallenged.txt.gz', sep = '')
    if (length(list.files(dir.data, pattern = nfit)) == 1){
      nfit <- paste(dir.data, nfit, sep = '')
      D <- read.csv(nfit, sep = ' ')
      return(D)
    } else if (length(list.files(dir.data, pattern = nfit)) == 0){
      print(paste(cllx, timepointx, 'no DEG list found', sep = ' '))
      D <- NULL
      return(D)
    } else if(length(list.files(dir.data, pattern = nfit)) >= 2){
      print(paste(cllx, timepointx, 'multiple DEG lists found?', sep = ' '))
      D <- NULL
      return(D)
    }
  }
  
  getAdjustedMean <- function(y) {
    m=zeroinfl(y~1|1, dist="negbin", link="logit", method ="Nelder-Mead") 
    mu=exp(m$coefficients$count)
    return(mu)
  }
  
  f_fc <- function(x, z1, z2, i) {
    zz1 <- as.numeric(x[i,which(colnames(x) %in% z1)])
    zz2 <- as.numeric(x[i,which(colnames(x) %in% z2)])
    if (min(zz1)==0 &
        min(zz2)==0) {
      if (max(zz1)!=0 &
          max(zz2)!=0){
        mean1 <- getAdjustedMean(zz1)
        mean2 <- getAdjustedMean(zz2)
      } else if (max(zz1)==0 &
                 max(zz2)!=0){
        mean1 <- 0
        mean2 <- getAdjustedMean(zz2)
      } else if (max(zz1)!=0 &
                 max(zz2)==0){
        mean1 <- getAdjustedMean(zz1)
        mean2 <- 0
      }
    } else if (min(zz1)!=0 &
               min(zz2)==0) {
      mean1 <- mean(zz1)
      mean2 <- getAdjustedMean(zz2)
    } else if (min(zz1)==0 &
               min(zz2)!=0) {
      mean1 <- getAdjustedMean(zz1)
      mean2 <- mean(zz2)
    } else if (min(zz1)!=0 &
               min(zz2)!=0) {
      mean1 <- mean(zz1)
      mean2 <- mean(zz2)
    }
    mean1 <- as.numeric(mean1)
    mean2 <- as.numeric(mean2)
    fc <- as.numeric(log10(mean1)-
      log10(mean2))
    return(list(fc, mean1, mean2))
  }
  
  f_enz2sym <- function(x){
    ref <- read.csv('data/gene_info_taxid-9606_version-20170511.tsv.gz', 
                    sep = '\t')
    ref <- ref[,2:3]
    # head(ref)
    gl <- ref[ref$GeneID %in% x,]
    return(gl)
  }
  
  f_sort <- function(x){
    gl <- x
    if (any(is.na(x$Symbol)) == TRUE){
      gl <- gl[!is.na(gl$Symbol),]
      print("NAs removed")
    } else {
      print("There are no NAs")
    }
    if (any(duplicated(gl$Symbol)) == TRUE){
      print(paste(length(which(duplicated(gl$Symbol))), "Duplicates removed", sep = ' '))
      gl <- gl[!duplicated(gl$Symbol),]
    } else {
      print("There are no duplicates")
    }
    return(gl)
  }
  
  cores <- detectCores()
  cl <- as.integer(cores[1]*0.9) # 90% of cores used
  registerDoParallel(cl)
  
  ######### Run FC calculations #########
  for (cll in 1:length(clls)){
    cllx <- clls[cll]
    print(cllx)
    for (timepoint in 1:length(timepoints)){
      timepointx <- timepoints[timepoint]
      print(timepointx)
      # Define the groups
      z <- colnames(X[,c(grep(cllx, colnames(X)))])
      z <- z[c(grep(timepointx, z))]
      z1 <- z[c(grep('_A_', z))] # Update dependent on group identifier
      z2 <- z[c(grep('_D_', z))] # Update dependent on group identifier
      # Load the DEGs
      D <- f_deg(cllx, timepointx)
      if (is.null(D)) { 
        print(cllx)
        next
      }
      # subset significant DEGs
      D <- D[D$qval < 0.05,]
      if (length(rownames(D))==0) { 
        print(paste(cllx, timepointx, 'has no significant DEGs', sep = ' '))
        next
      }
      reps <- intersect(D$gene_short_name, rownames(X))
      Y <- X[which(rownames(X) %in% reps),]
      Y <- Y[,which(colnames(Y) %in% z)]
      # Calculate the FCs, log10(mean(z1))-log10(mean(z2))
      fc_val <- foreach (i=1:length(D$gene_short_name), .combine = cbind) %dopar% {
        tempd <- f_fc(Y, z1, z2, i)
      }
      # Add FC to DEG files. Note, defined colnames are analysis dependent
      if(length(colnames(fc_val)) >= 1){
        D$'FC' <- fc_val[1,]
        D$'mean_expr_challenged' <- fc_val[2,] # Update dependent on group identifier
        D$'mean_expr_unchallenged' <- fc_val[3,] # Update dependent on group identifier
      } else if (length(colnames(fc_val)) == 0) {
        D$'FC' <- fc_val[1]
        D$'mean_expr_challenged' <- fc_val[2] # Update dependent on group identifier
        D$'mean_expr_unchallenged' <- fc_val[3] # Update dependent on group identifier
      }
      ref <- f_enz2sym(D$gene_short_name) 
      D <- left_join(D, ref, by = c('gene_short_name' = 'GeneID'))
      D <- f_sort(D) ## Duplicates and zeros are here removed
      nFIT <- paste(dir.out, 
                    "/degs_monocle_", type,  '_', cllx, '_', timepointx, 
                    "_AllergenChallenged_vs_NonChallenged.txt", sep = '') # Update dependent on group identifier
      outf <- D[,c('Symbol', 'FC', 'qval', 'mean_expr_challenged', 'mean_expr_unchallenged')] # Update dependent on group identifier
      outf$FC <- as.numeric(outf$FC)
      outf$mean_expr_challenged <- as.numeric(outf$mean_expr_challenged) # Update dependent on group identifier
      outf$mean_expr_unchallenged <- as.numeric(outf$mean_expr_unchallenged) # Update dependent on group identifier
      write.table(outf, nFIT, sep = '\t', row.names = F, col.names = T, quote = F)
    }
  }
}
