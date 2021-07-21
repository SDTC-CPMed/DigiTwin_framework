# scRNA-seq
## Extraction of sequencing output data

1. getSeq_3samples.sh # Data extraction
2. getDgem.sh # Create digital gene expression matrix and fastqc analysis
3. getdepth.sh # Check read depth

## Quality assesment and full matrix construction

1. sc_data_quality_sorting.R # Check the quality of the data and remove outliers
2. full_sc-matrix_construction.R # Combine the matrices from all samples and set specific quality cut-offs

## Cell type analysis

1. RCA_reference_construction.R # Build references for cell type identification
2. RCA_cellType_identification.R # Cell type identification using reference component analysis

## Differnetial expression analysis 

1. pre-DEG_analysis.R # statistics and preparation of input files for DEG analysis  
2. Monocle_v3_RAdata.r # Calculate DEGs
3. sc_FC_zero-infl-neg-bionmial.R # Fold change calculations
4. plot_DEGs.R # Plot the distribution of DEGs
