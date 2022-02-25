#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# set default parameters
time_points <- c('0h', '12h', '1D', '2D', '3D', '5D', '7D')

stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

if (args[1]=='--help') {
  cat('singleR.R <full path to ipa-input directory> <full path to deg-input directory> <full output path>\n')
  cat('Input ipa directory:
        Should contain subdirectories, one for each time point (time point need to be in subdirectory name). 
        Ensure that the directory  only contains files included in this analysis, or error will occur when listing files.
        Ensure the names of files follows the system <[cellType]_[timePoint]_[other description].xls>\n')
  cat('Input deg directory:
        Ensure time points and cell type are named equally as in the ipa data file names.
        The files should be named <degs_monocle_HA_vs_HC_[CellType]_[timePoint].txt> (Note, CellType in 6th position sep="_")\n')
  cat('Output directory: Will be created if it does not exist.
      If output files already exists, a number will be added to its name to prevent overwriting\n')
  cat('optional parameters;\n')
  cat('--time_points="str1" "str2" ... (default: 0h 12h 1D 2D 3D 5D 7D)\n') 
  cat('Built on R version 4.1.2\n')
  stop_quietly()
}

# Test so that the number of arguments are correct: if not, return an error
if (length(args)<3) {
  stop("At least three arguments must be supplied (see --help guide).n", call.=FALSE)
}
# read optional parameters
if (length(args)>3) {
  if (length(grep('--time_points', args, value = TRUE))==1){
    analysis = 'time_points'
    time_points <- c(strsplit(grep('--time_points', args, value = TRUE), split = '=')[[1]][[2]])
  } else {
    stop(paste("Don't recognize some of the flags;", args[3:length(args)]), call.=FALSE)
  }
}


library(readxl)
library(dplyr)

# dir.home <- getwd()
# dir.data <- paste(dir.home, '/example_data/IPA_UR-prediction', sep = '')
# dir.degs <- paste(dir.home, '/example_data/DEGs', sep = '')
# dir.out <- paste(dir.home, '/output/test', sep = '')
dir.data <- args[1]
dir.degs <- args[2]
dir.out <- args[3]

if (dir.exists(dir.out)==FALSE){
  dir.create(dir.out, recursive = T)
  print(paste('created', dir.out))
}

molecule_types <- c('complex', 'cytokine', 'G-protein coupled receptor',
                    'group', 'growth factor', 'ligand-dependent nuclear receptor',
                    'transmembrane receptor') 
ipa_columns <- c("Upstream Regulator", "Expr Log Ratio", "Molecule Type",
                 "Predicted Activation State", "Activation z-score",
                 "p-value of overlap", "Target Molecules in Dataset",
                 "Mechanistic Network")

# for each time point, create the MNM
for (time_point in time_points){
  print(time_point)
  # time_point <- time_points[5]
  # read in the ipa output files
  print('Load the IPA data')
  indir <- list.files(dir.data, pattern = time_point, full.names = T)
  ipa_files <- list.files(indir, full.names = T)
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
  celltypes <- list.files(indir)
  celltypes <- gsub(' ', '', celltypes)
  celltypes <- sapply(strsplit(celltypes, '_'), '[[', 1)
  names(ipa) <- celltypes
  
  # read in the DEG files
  print('Load the DEG data')
  deg_files <- list.files(dir.degs, pattern = time_point, full.names = T) 
  degs <- c()
  for (i in 1:length(deg_files)){
    degs[[i]] <- read.table(deg_files[i], header = T)
    # Only include positive FCs
    degs[[i]] <- degs[[i]][which(degs[[i]]$FC > 0),]
  }
  celltypes <- list.files(dir.degs, pattern = time_point) 
  celltypes <- sapply(strsplit(celltypes, '_'), '[[', 6)
  names(degs) <- celltypes
  
  # Create MNM
  print('Create MNM')
  for (i in 1:length(ipa)){
    # i <- 1
    target_ct <- names(ipa)[i]
    urs <- ipa[[i]]$`Upstream Regulator`
    for (z in 1:length(degs)){
      # z <- 1
      source_ct <- names(degs)[z]
      # check which urs are among deg
      urs_deg <- urs[which(urs %in% degs[[z]]$Symbol)]
      if (any(grep('PDGF BB', urs))){
        if ('PDGFB' %in% degs[[z]]$Symbol){
          urs_deg <- c(urs_deg, 'PDGFB')
        }
      }
      # construct output file
      ipa_sub <- ipa[[i]][which(ipa[[i]]$`Upstream Regulator` %in% urs_deg),]
      degs_sub <- degs[[z]][which(degs[[z]]$Symbol %in% urs_deg),]
      ipa_sub <- full_join(ipa_sub, degs_sub, by = c('Upstream Regulator' = 'Symbol'))
      ipa_sub$'Source cell type' <- source_ct
      ipa_sub$'Target cell type' <- target_ct
      if (i == 1 & z == 1){
        MNM_out <- ipa_sub
      } else {
        colnames(MNM_out) == colnames(ipa_sub)
        MNM_out <- rbind(MNM_out, ipa_sub)
      }
    }
  }
  
  # Save output file
  print('Save to output')
  out_name <- paste(time_point, '_UR_interactions.csv', sep = '')
  # Control not to overwrite existing files
  if (file.exists(paste(dir.out, out_name, sep = '/'))){
    out_name <- paste(sapply(strsplit(out_name, '.csv'), '[[', 1), '_',
                      length(list.files(dir.out, pattern = out_name)), '.csv', sep = '')
    print('Output file already exists. Name adjusted not to overwrite')
  }
  write.csv(MNM_out, paste(dir.out, out_name, sep = '/'), row.names = F)
  print(paste(out_name, 'was created'))
  
}

