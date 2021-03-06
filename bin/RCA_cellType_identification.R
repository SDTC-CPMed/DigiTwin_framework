# Sandra Lilja

#' Cell type identification using reference component analysis
#' 
#' @param tissue Either 'HA' for allergic patients of 'HC' for healthy controls
#' 
#' @return no return, just saves the files into data folder
#' @export
#'


library(dplyr)
library(RCA)
library(data.table)
library(edgeR)
library(R.utils)
library(picante)


RCA_cellType_identification = function(tissue){
  if(tissue == 'HC'){k = 'k21'
  } else{k = 'k14'}
  
  dir.data <- 'data/knn_smoothing/'
  dir.out <- 'data/RCA_out/'
  
  
  featureConstructMy <- dget("featureConstruct_pvals.R")
  cellClustMy <- dget("cellClust.R")
  RCAPlotMy <- dget("RCAPlot.R")
  get_RCA_labels_and_PVals <- dget("get_RCA_labels_and_PVals.r")
  
  data(sysdata, envir = environment())
  names(sysdata)

  
  for (roundx in 1:3){
    if (roundx == 1){
      ### Round 1
      X <- read.csv(paste(dir.data, tissue, '_min200genesPerCell_sorted_expression_matrix.knn-smooth_', k, '.csv.gz', sep = ''), sep = ',', header = T, row.names = 1)

    } else if (roundx == 2){
      ### Round 2
      clust.method <- 'GlobalPanel_MonvsDCvsBvsTNK'
      X <- readRDS(paste(dir.out, tissue, '_min200genesPerCell_sorted_expression_matrix.knn-smooth_', k, '_und-rm_withCellTypes_', clust.method, '.rds', sep = ''))
      
      X[1:5,1:5]
      X <- X[,grep('TNK_', colnames(X))]
      colnames(X) <- sapply(strsplit(colnames(X), 'NK_'), '[[', 2)
    } else if (roundx == 3){
      # ### Round 3
      clust.method <- 'CellRef_Huan_PBMC_CD4vsCD8vsNK'
      X <- readRDS(paste(dir.out, tissue, '_min200genesPerCell_sorted_expression_matrix.knn-smooth_', k, '_und-rm_withCellTypes_', clust.method, '.rds', sep = ''))
      
      X[1:5,1:5]
      X <- X[,grep('CD4_', colnames(X))]
      colnames(X) <- sapply(strsplit(colnames(X), 'D4_'), '[[', 2)
    }
    ########################################################
    X[1:5,1:5]
    ############## pre RCA ##############
    ### Define the input data set
    X <- as.data.frame(X)
    
    f.rm.und <- function(A, D, p = 0.05){
      get.RCA.labels <- dget("get_RCA_labels_and_PVals.r")
      labels = get.RCA.labels(D)
      z <- (is.na(labels$RCA_Pval_adj)) | (labels$RCA_Pval_adj >= p)
      und <- rownames(labels[z,])
      A <- A[, !(colnames(A) %in% und)]
      list(expr=A, und=und)
    }
    

    ENTREZ.origin <- X


    Xin <- X
    X[1:5,1:5]
    
    #######################################
    #### Convert GeneSymbols to entrez ####
    #######################################
    ref <- read.csv('data/gene_info_taxid-9606_version-20170511.tsv.gz', 
                    sep = '\t')
    ref <- ref[,2:3]
    head(ref)
    
    Xin[,'Symbol'] <- rownames(Xin)
    Xin[1:5,(length(colnames(Xin))-5):length(colnames(Xin))]
    Xin <- left_join(Xin, ref)
    Xin[1:5,(length(colnames(Xin))-5):length(colnames(Xin))]
    Xin <- Xin[!is.na(Xin$GeneID),]
    rownames(Xin) <- Xin$GeneID
    Xin <- Xin[,-length(colnames(Xin))]
    Xin <- Xin[,-length(colnames(Xin))]
    Xin[1:5,1:5]
    
    remove(ref)
    #################################################
    
    f.data.obj <- function(Y, case, logged, clust.method=clust.method){
      data.obj <- dataConstruct(Y)
      data.obj <- geneFilt(obj_in = data.obj)
      data.obj <- cellNormalize(data.obj) # default "no_norm"
      data.obj <- dataTransform(data.obj, 'log10') # change values (fpkm_transformed)
      data.obj <- featureConstructMy(data.obj, method = clust.method) # Change values (fpkm_for_clust)
      data.obj <- cellClustMy(data.obj, deepSplit_wgcna = 1, min_group_Size_wgcna = 5)
    }
    
    get.RCA.labels <- dget("get_RCA_labels_and_PVals.r")
    
    names(sysdata)
    und.rm <- T
    X <- Xin
    
    # RCA analysis
    if (roundx == 1){
      clust.method <- "GlobalPanel_MonvsDCvsBvsTNK"
    } else if (roundx == 2){
      clust.method <- "CellRef_Huan_PBMC_CD4vsCD8vsNK"
    } else if (roundx == 3){
      clust.method <- "CellRef_Huan_PBMC_CD4subsets"
    }
    
    Xin <- Xin[intersect(rownames(as.data.frame(sysdata[clust.method])), rownames(Xin)),]
    X <- X[,colSums(Xin)!=0]
    ###############################################
    ### Change all 1 to 1.1 in dataset to solve issue with log(1)=0 after transformation
    X[X==1] <- 1.1
    X <- as.data.frame(X)
    ## run RCA
    data.obj <- dataConstruct(X)
    data.obj <- geneFilt(obj_in = data.obj)
    data.obj <- cellNormalize(data.obj) # default "no_norm"
    data.obj <- dataTransform(data.obj, 'log10') # change values (fpkm_transformed)
    data.obj <- featureConstructMy(data.obj, method = clust.method) # Change values (fpkm_for_clust)
    data.obj <- cellClustMy(data.obj, deepSplit_wgcna = 1, min_group_Size_wgcna = 5)
    Y <- X
    # subset cells which match their respective reference cell type. 
    data.obj_backup = data.obj
    
    if(und.rm){
      und <- data.obj$und
      i <- 0
      Y <- data.obj$fpkm_raw
      repeat{
        aux <- f.rm.und(Y, data.obj)
        n <- length(aux[["und"]])
        if(n == 0){
          break
        }
        Y <- aux[["expr"]]
        und <- c(und, aux[["und"]])
        data.obj <- f.data.obj(Y, case, logged, clust.method=clust.method)
        und <- c(und, data.obj$und)
        i <- i + 1
        cat("it-", i, " und: ", length(und), "\n")
        cat(dim(Y))
        cat("\n")
      }
    }
    
    # extract data for t-sne plots
    fpkm_for_clust <- data.obj$fpkm_for_clust
    corrScorRHO <- data.obj$corrScorRHO
    
    labels = get.RCA.labels(data.obj)
    labels <- labels[intersect(colnames(ENTREZ.origin), rownames(labels)),]
    samplename <- tissue
    write.table(labels, paste(dir.out, samplename, '_RCA_cell_types_', clust.method, '.txt', sep = ''), sep = '\t', col.names = T, row.names = T, quote = F)
    write.csv(fpkm_for_clust, paste(dir.out, samplename, '_fpkm_for_clust_', clust.method, '.csv', sep = ''))
    write.csv(corrScorRHO, paste(dir.out, samplename, '_corrScorRHO_', clust.method, '.csv', sep = ''))
    

    RCAPlotMy(data.obj)
    files <- list.files(dir.out, pattern = '^RCAplot*', full.names = T)
    samplename <- paste(tissue, "und-rm", clust.method, sep = '_')
    sapply(files, FUN = function(eachPath){
      file.rename(from = eachPath, to= sub(pattern = paste('\\', dir.out, '/', sep = ''), paste0(paste('\\', dir.out, '/', samplename, '_', sep = '')), eachPath))
    })
    
    ### Rename and sort out cells in the original matrix according to new IDs
    out.X <- ENTREZ.origin[,rownames(labels)]
    colnames(out.X) <- paste(gsub('_', '-', labels$RCA), colnames(out.X), sep = '_')
    length(grep(FALSE, gsub('-', '_', sapply(strsplit(colnames(out.X), '_'), '[', 1)) == labels$RCA))
    saveRDS(out.X, file = paste(dir.out, tissue, '_min200genesPerCell_sorted_expression_matrix.knn-smooth_', k, '_und-rm_withCellTypes_', clust.method, '.rds', sep = ''))
    
    
    # Combine all the results into one matrix
    if (roundx == 3){
      ### for PBMC
      ### create full matrix
      names(sysdata)
      clust.method <- "GlobalPanel_MonvsDCvsBvsTNK"
      X1 <- readRDS(paste(dir.out, tissue, '_min200genesPerCell_sorted_expression_matrix.knn-smooth_', k, '_und-rm_withCellTypes_', clust.method, '.rds', sep = ''))
      
      clust.method <- "CellRef_Huan_PBMC_CD4vsCD8vsNK"
      X2 <- readRDS(paste(dir.out, tissue, '_min200genesPerCell_sorted_expression_matrix.knn-smooth_', k, '_und-rm_withCellTypes_', clust.method, '.rds', sep = ''))
      
      clust.method <- "CellRef_Huan_PBMC_CD4subsets"
      X3 <- readRDS(paste(dir.out, tissue, '_min200genesPerCell_sorted_expression_matrix.knn-smooth_', k, '_und-rm_withCellTypes_', clust.method, '.rds', sep = ''))
      
      
      # remove sub-clustered cell types
      X1[1:5,1:5]
      X1 <- X1[,!grepl('TNK_', colnames(X1))]
      X2[1:5,1:5]
      X2 <- X2[,!grepl('CD4_', colnames(X2))] 
      X3[1:5,1:5]
      
      # combine the data
      X <- cbind(X1,X2,X3)
      X[1:5,1:5]
      
      ref <- read.csv('data/gene_info_taxid-9606_version-20170511.tsv.gz', sep = '\t')
      
      ref <- ref[,2:3]
      head(ref)
      
      X[,'Symbol'] <- rownames(X)
      X[1:5,(length(colnames(X))-5):length(colnames(X))]
      X <- left_join(X, ref)
      X[1:5,(length(colnames(X))-5):length(colnames(X))]
      X <- X[!is.na(X$GeneID),]
      rownames(X) <- X$GeneID
      X <- X[,-length(colnames(X))]
      X <- X[,-length(colnames(X))]
      X[1:5,1:5]
      

      saveRDS(X, file = paste(dir.out, 'full_matrix/', tissue, '_min200genesPerCell_sorted_ENTREZ_expression_matrix.knn-smooth_', k, '_und-rm_withCellTypes.rds', sep = ''))
      
    }
  }
  
} 

