Below scripts works using R 3.4 and python 3.7.4.  

Note that the R scripts are not developed to be executed from command line. We recommend using RStudio for running these scripts. Modifications to these codes may be needed, to set input file directories, sample names, and other project specific information.  

# scRNA-seq
## Quality assesment and full matrix construction

To ensure good quality data for downstream analyses, it is recommended to remove poor quality cells from the analyses. Here, we define and keep the good quality cells as those having a minimum of 400 transcripts, 200 genes, and less than 20% mitochondrial genes. Due to the risk of duplicates in the library resulting in two or more cells sharing a cell barcode, it is also recommend to remove outliers, which can be based on empirical evaluation of the distributionan overestimation of transcripts count over the cells. 

First, we combined the output files from getDgem into a matrix using full_sc-matrix_construction.R including some minimum quality cut-offs, whereafter we estimate the overall quality of the data and remove outliers using sc_data_quality_sorting.R. 

1. full_sc-matrix_construction.R # Combine the matrices from all samples and set specific quality cut-offs
2. sc_data_quality_sorting.R # Check the quality of the data and remove outliers

## Cell type analysis

The references for cell type identification were created by RCA_reference_construction.R in R 3.4, and added to the RCA sysdata file installed. Input to this scripts are the normalized matrix output from microarray analyses based on which the references should be constructed.   

RCA_cellType_identification.R is then run to map the cells towards the newly created references. As this script is based on some modifications to the Reference Component Analysis (RCA, Li et al., 2017 https://pubmed.ncbi.nlm.nih.gov/28319088/) with calculations of p-values and newly created references, the modified codes need to be added as well (as part of the script). These modified codes you can find in /RNA_mod/. As input matrix to this code we used the knn-smoothed data (output using script from Wagner et al., 2018 https://www.biorxiv.org/content/early/2018/04/09/217737), but any UMI-based expression matrix should work equally.  

1. RCA_reference_construction.R # Build references for cell type identification
2. RCA_cellType_identification.R # Cell type identification using reference component analysis

## Differnetial expression analysis 

To prepare the data from cell type analysis for differential expression analysis, we run pre-DEG_analysis.R. This script will ensure that the input to Monocle_v3_RAdata.r is is the correct format. It will also compute some basic statistics, such as the number of cells per sample (cell type, time point, etc). 

Monocle_v3_RAdata.r is run to calculate differentially expressed genes (DEGs) between allergen challenged and non-challenged samples. Some modifications to the code are needed to compute DEGs between any other pair of groups.

The fold changes (FCs) for the DEGs, we calculate using sc_FC_zero-infl-neg-bionmial.R. This will compute FCs between allergen challenged vs non-challenged, meaning that a positive FC indicates a higher mean expression in the allergen challenged compared to the non-challenged group. Some modifications to the code are needed to compute the FCs between any other pair of groups.

We plot the distrinbution of DEGs over different time points and cell types (Fig 1), as well as the number of cell types over different time points and treatment groups (Fig 2), by plot_DEGs.R. 

1. pre-DEG_analysis.R # statistics and preparation of input files for DEG analysis  
2. Monocle_v3_RAdata.r # Calculate DEGs
3. sc_FC_zero-infl-neg-bionmial.R # FC calculations
4. plot_DEGs.R # Plot the distribution of DEGs

<img src="https://user-images.githubusercontent.com/51739216/144604687-54ee4a7c-b661-4ec7-927b-49d54efe2883.png" width="500" />  
Fig 1. log(number of DEGs) identified between allergen-stimulated and diluent-stimulated cells, for each cell type at the different time points.  

<img src="https://user-images.githubusercontent.com/51739216/144604759-9677b860-46a9-4296-876b-eb5de3b69b1a.png" width="500" />  
Fig 2. Cell type proportions in different groups of allergen-stimulated (A) and diluent-stimulated (D), as well as unstimulated control (C), samples from healthy individuals at the different time points.  
