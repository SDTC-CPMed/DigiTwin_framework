# Sandra Lilja

library(dplyr)
library(plyr)
library(ggplot2)
library(scales)
library(reshape2)


plot_custom = function(){

  dir.data <- 'data/DEGS_with_Monocle/Monocle_out_withFCs/'
  dir.out.images <- 'plot/'
  
  #### Count all DEGs for each cell type and timepoint
  nm <- list.files(path=dir.data)
  cells <- unique(sapply(strsplit(list.files(path=dir.data), '_'), '[[', 4))
  timepoint <- unique(sapply(strsplit(list.files(path=dir.data), '_'), '[[', 5))
  X <- matrix(data = NA, nrow = length(cells), ncol = length(timepoint))
  colnames(X) <- timepoint
  rownames(X) <- cells
  
  for (i in 1:length(nm)){
    Y <- read.table(paste(dir.data, nm[i], sep = ''), header = T, sep = '\t')
    nDEGs <-  length(Y[,1])
    a <- unique(sapply(strsplit(nm[i], '_'), '[[', 4))
    b <- unique(sapply(strsplit(nm[i], '_'), '[[', 5))
    X[a,b] <- nDEGs
    remove(Y, nDEGs, a, b)
  }
  X[is.na(X)] <- 0
  
  #### Plot nDEGs over time for each cell type
  colnames(X) <- c('day0.5', 'day1', 'day2', 'day3', 'day5', 'day7')
  X <- as.data.frame(X)
  X <- log(X)
  X[X=='-Inf'] <- 0
  X$celltype <- as.character(rownames(X))
  ndata <- melt(X, id = "celltype")
  ndata <- arrange(ndata, celltype)
  ndata <- mutate(ndata, variable = gsub("day", "", ndata$variable))
  colnames(ndata) = c('celltype', 'day', 'log_nDEGs')
  #ndata <- rename(ndata, day = variable)
  #ndata <- rename(ndata, log_nDEGs = value)
  ndata <- mutate(ndata, day = as.numeric(day))
  newplot <- ggplot(ndata, aes(x = day, y = log_nDEGs, color = celltype, group = celltype)) +
    geom_point() +
    geom_line()
  pdf(paste(dir.out.images, '/lognDEGs_over_time3.pdf', sep = ''))
  newplot
  dev.off()
  
  ##############################################################################################
  #### Plot the distribution of different cell types for each separate timepoint and treatment
  library(plyr)
  library(ggplot2)
  library(scales)
  library(reshape2)
  #rm(list=ls())
  
  #dir.home <- getwd()
  dir.data <- 'data/RCA_out/full_matrix/'
  
#  aa <- list.files(dir.data, pattern = 'ENTREZ')
#  aa <- list.files(dir.data)
  
  #X <- read.table(paste(dir.data, aa, sep = ''), header = T, sep = '\t', row.names = 1)
  X <- readRDS('data/RCA_out/full_matrix/HA_min200genesPerCell_sorted_ENTREZ_expression_matrix.knn-smooth_k14_und-rm_withCellTypes.rds')
  
  X[1:5,1:5]
  
  M <- matrix(NA, ncol = 3, nrow = length(colnames(X)))
  M[,1] <- sapply(strsplit(colnames(X), '_'), '[[', 1)
  M[,2] <- sapply(strsplit(colnames(X), '_'), '[[', 3)
  M[,3] <- sapply(strsplit(colnames(X), '_'), '[[', 4)
  colnames(M) <- c('CellType', 'timepoint', 'treatment')
  head(M)
  ### Create data frame 'Mdf' with unique rows containing column 'value' 
  ### specifying the ratio of the cell type in its certain timepoint and condition.
  df <- as.data.frame(M)
  df <- ddply(df,.(CellType, timepoint, treatment),nrow)
  colnames(df) <- c('CellType', 'timepoint', 'treatment', 'ratio')
  
  newplot <- ggplot(df, aes(x=treatment, y=ratio, fill=CellType)) +
    geom_bar(stat="identity", position = "fill") +
    scale_y_continuous(labels = percent_format())+
    facet_grid(~timepoint) + theme_bw()
  pdf(paste(dir.out.images, '/celltype_ratios_over_groups4.pdf', sep = ''))
  
  newplot
  dev.off()
  
  
}

