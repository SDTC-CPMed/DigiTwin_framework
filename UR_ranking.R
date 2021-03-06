#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

if (args[1]=='--help') {
  print('singleR.R <full path to ipa-input main directory (i.e., containing subdirectories and data from all samples to be compared)> <full output path>')
  print('Built on R version 4.1.2')
  stop_quietly()
}

# Test so that the number of arguments are correct: if not, return an error
if (length(args)<2) {
  stop("At least two arguments must be supplied (see --help guide).n", call.=FALSE)
}


library(readxl)
library(dplyr)
library(pheatmap)

# dir.home <- getwd()
# dir.data <- paste(dir.home, '/ipa_output', sep = '')
# dir.out <- paste(dir.home, '/output', sep = '')
dir.data <- args[1]
dir.out <- args[3]

if (dir.exists(dir.out)==FALSE){
  dir.create(dir.out, recursive = T)
  print(paste('created', dir.out))
}

molecule_types <- c('complex', 'cytokine', 'G-protein coupled receptor',
                    'group', 'growth factor', 'ligand-dependent nuclear receptor',
                    'transmembrane receptor') 
ipa_columns <- c("Upstream Regulator", "Molecule Type",
                 "Activation z-score", "p-value of overlap")


# Rank the URs based on in how many cell types and at how many time points they were predicted in

# read in the ipa output files
ipa_files <- list.files(dir.data, recursive = T, pattern = '^[^UR_]', full.names = T)
ipa_files <- ipa_files[-grep('manually_processed', ipa_files)]

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
celltypes <- sapply(strsplit(ipa_files, '/'), '[[', 9)
celltypes <- gsub(' ', '', celltypes)
time_points <- sapply(strsplit(celltypes, '_'), '[[', 2)
celltypes <- sapply(strsplit(celltypes, '_'), '[[', 1)
celltypes <- gsub('monocytes', 'Monocytes', celltypes)
celltypes <- gsub('th1', 'Th1', celltypes)
celltypes <- gsub('Bcells', 'Bcell', celltypes)
celltypes <- gsub('Th1cell', 'Th1', celltypes)
celltypes <- gsub('Th17cell', 'Th17', celltypes)
celltypes <- gsub('Th2cell', 'Th2', celltypes)
celltypes <- gsub('CD8cell', 'CD8', celltypes)
celltypes <- gsub('Tregcell', 'Treg', celltypes)
celltypes <- gsub('Monocytescell', 'Monocytes', celltypes)
celltypes <- gsub('Dendriticcell', 'Dendritic', celltypes)
celltypes <- gsub('NKcell', 'NK', celltypes)
celltypes <- gsub('NTcell', 'NT', celltypes)
celltypes <- gsub('Bcell', 'Bcells', celltypes)
unique(sort(celltypes))

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

write.csv(ipa_comb, paste(dir.out, '/combined_UR_regulations.csv', sep = ''), row.names = F)
write.csv(UR_rank, paste(dir.out, '/UR_ranking.csv', sep = ''), row.names = F)



# Make a heatmap of the z-scores for each UR in each cell type at each time point
zscores <- matrix(NA, ncol = length(ipa), nrow = length(URs))
colnames(zscores) <- names(ipa)
rownames(zscores) <- UR_rank$UR
for (i in 1:length(colnames(zscores))){
  for (z in 1:length(rownames(zscores))){
    # i <- 1
    # z <- 1
    subx <- ipa_comb[which(ipa_comb$Sample == colnames(zscores)[i]),]
    subx <- subx[which(subx$`Upstream Regulator` == rownames(zscores)[z]),]
    if (length(subx$`Activation z-score`) == 1){
      zscores[z,i] <- subx$`Activation z-score`
    } else if (length(subx$`Activation z-score`) == 0){
      zscores[z,i] <- 0
    } else if (length(subx$`Activation z-score`) > 1){
      print(paste(i, z, 'Error: more than one z-score selected'))
    }
  }
}

# create the heatmap for FC data
bk1 <- c(seq(-10,-0.1,by=0.1),-0.001)
bk2 <- c(0.001,seq(0.1,10,by=0.1))
bk <- c(bk1,bk2)  #combine the break limits for purpose of graphing

my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(bk1)-1),
                "black", "black",
                c(colorRampPalette(colors = c("pink", "darkred"))(n = length(bk2)-1)))

annotation_colors <- data.frame(row.names = colnames(zscores), 
                                CellType = sapply(strsplit(colnames(zscores), '_'), '[[', 1),
                                TimePoint = sapply(strsplit(colnames(zscores), '_'), '[[', 2))

ph <- pheatmap(zscores, color = my_palette, breaks = bk,
         show_colnames = F,
         annotation_col = annotation_colors, 
         cluster_cols = F,
         cluster_rows = F,
         cellwidth = 4)

pdf(paste(dir.out, '/UR_ranking_heatmap.pdf', sep = ''), height = 13, width = 8)
print(ph)
dev.off()
