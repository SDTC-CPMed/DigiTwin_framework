Below scripts works using R 3.4 and python 3.7.4. 

# scRNA-seq
## Extraction of sequencing output data

1. getSeq_3samples.sh # Data extraction
2. getDgem.sh # Create digital gene expression matrix and fastqc analysis
3. getdepth.sh # Check read depth
4. example_data_extraction.sh # example script to run above scripts

Run examples:  
getSeq_3samples [path/input_file] [path_to_GenomeDir] [path/primary_assembly.genome.fa] [path/primary_assembly.annotation.gtf] 20000  
where 20000 is the double number of cells expected to be extracted from analysis.

getDgem [path/output_directory] [path_input_data] [run_ID] [sample_ID]  
getdepth [path/output_directory] [path_input_data] [run_ID] [sample_ID]  
where sample_ID has to be the same as output sample name from getSeq_3samples

See example_data_extraction.sh for an example on how to extract data from one sequencing run. 

## Quality assesment and full matrix construction

To ensure good quality data for downstream analyses, it is recommended to remove poor quality cells from the analyses. Here, we define and keep the good quality cells as those having a minimum of 400 transcripts, 200 genes, and less than 20% mitochondrial genes. Due to the risk of duplicates in the library resulting in two or more cells sharing a cell barcode, it is also recommend to remove outliers, which can be based on empirical evaluation of the distributionan overestimation of transcripts count over the cells. 

First, we combined the output files from getDgem into a matrix using full_sc-matrix_construction.R including some minimum quality cut-offs, whereafter we estimate the overall quality of the data and remove outliers using sc_data_quality_sorting.R. 

1. full_sc-matrix_construction.R # Combine the matrices from all samples and set specific quality cut-offs
2. sc_data_quality_sorting.R # Check the quality of the data and remove outliers

## Cell type analysis

1. RCA_reference_construction.R # Build references for cell type identification
2. RCA_cellType_identification.R # Cell type identification using reference component analysis

## Differnetial expression analysis 

1. pre-DEG_analysis.R # statistics and preparation of input files for DEG analysis  
2. Monocle_v3_RAdata.r # Calculate DEGs
3. sc_FC_zero-infl-neg-bionmial.R # Fold change calculations
4. plot_DEGs.R # Plot the distribution of DEGs
