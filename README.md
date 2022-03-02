README
================

Below scripts works using R 3.4, R 4.0 and python 3.7.4.

Note that the R scripts are not developed to be executed from command
line. We recommend using RStudio for running these scripts.
Modifications to these codes may be needed, to set input file
directories, sample names, and other project specific information.

# scRNA-seq analysis for construction of multicellular network models (MNM) and prioritization of upstream regulatory (UR) genes

# Instructions for running the codes

The project is divided into 3 main parts

-   [Data preparation, and analyses for UR prediction and MNM
    construction](#data-preparation-and-analyses-for-ur-prediction-and-mnm-construction),
    which can be run by the script bin/main.R.

-   [Ingenuity Pathway
    analysis](#upstream-regulatory-ur-gene-prediction-in-ingenuity-pathway-analysis-IPA).
    The detailed instructions for running the analyses follow below.

-   [Construction of multicellular network models (MNMs) and ranking of
    URs](#construction-of-multicellular-network-models-mnms-and-ranking-of-urs),
    which can be run by the script post_IPA.sh

## Data preparation, and analyses for UR prediction and MNM construction

The data can be processed by running the following script. The file
creates the directory structure and calls the functions to run [Quality
assesment and full matrix
construction](#quality-assesment-and-full-matrix-construction), [Cell
type analysis](#cell-type-analysis) and [Differnetial expression
analysis](#differnetial-expression-analysis)

``` eval
dir.create("data")
dir.create("data/DEGS_with_Monocle")
dir.create("data/DEGS_with_Monocle/Matrix_in")
dir.create("data/DEGS_with_Monocle/Monocle_out")
dir.create("data/DEGS_with_Monocle/Monocle_out_withFCs")
dir.create("data/knn_smoothing")
dir.create("data/RCA_out")
dir.create("data/RCA_out/full_matrix")
dir.create("data/RCA_references")
dir.create("plot")

### Read data from GSE180697 and save it to data folder
filenames = c('data/GSE180697_SAR_patients_expression_matrix.csv.gz',
              'data/GSE180697_Healthy_controls_expression_matrix.csv.gz')

### Preprocess the data and remove outliers
source('sc_data_quality_sorting.R')
HA = sc_data_quality_sorting(filenames[1])
write.csv(HA, 'data/HA_min200genesPerCell_sorted_expression_matrix.csv', quote = F)
HC = sc_data_quality_sorting(filenames[2])
write.csv(HC, 'data/HC_min200genesPerCell_sorted_expression_matrix.csv', quote = F)

### Run knn smoothing
#system(chmod u+x run_knn_smoothing.sh)
system('./run_knn_smoothing.sh ')

### Run cell typing
source('RCA_reference_construction.R')
source('RCA_cellType_identification.R')
RCA_cellType_identification('HA')
#RCA_cellType_identification('HC')

### Run pre-DEG analysis
source('pre-DEG_analysis.R')
pre_DEG_analysis('HA')
#pre_DEG_analysis('HC')

### Run Monocle_v3_RAdata
source('Monocle_v3_RAdata.r')
Monocle_v3_RAdata('HA')
#Monocle_v3_RAdata('HC')

### Run sc_FC_zero_infl_neg-binomial
source('sc_FC_zero-infl-neg-bionmial.R')
sc_FC_zero_infl_neg_binomial('HA')
#sc_FC_zero_infl_neg_binomial('HC')

### Run plot_DEGs.R
source('plot_DEGs.R')
plots = plot_custom()
pdf('plot/lognDEGs_over_time.pdf')
plots[[1]]
dev.off()
pdf('plot/celltype_ratios_over_groups.pdf')
plots[[2]]
dev.off()
```

### Quality assesment and full matrix construction

To ensure good quality data for downstream analyses, it is recommended
to remove poor quality cells from the analyses. Here, we define and keep
the good quality cells as those having a minimum of 400 transcripts, 200
genes, and less than 20% mitochondrial genes. Due to the risk of
duplicates in the library resulting in two or more cells sharing a cell
barcode, it is also recommend to remove outliers, which can be based on
empirical evaluation of the distributionan overestimation of transcripts
count over the cells.

The expected input is for this part is a matrix with genes in rows and
cells in columns. The quality of the data and removal of the outliers is
done using sc_data_quality_sorting.R

### Cell type analysis

The references for cell type identification were created by
RCA_reference_construction.R in R 3.4, and added to the RCA sysdata file
installed. Input to this scripts are the normalized matrix output from
microarray analyses based on which the references should be constructed.

RCA_cellType_identification.R is then run to map the cells towards the
newly created references. As this script is based on some modifications
to the Reference Component Analysis (RCA, Li et al., 2017
<https://pubmed.ncbi.nlm.nih.gov/28319088/>) with calculations of
p-values and newly created references, the modified codes need to be
added as well (as part of the script). These modified codes you can find
in /RNA_mod/. As input matrix to this code we used the knn-smoothed data
(output using script from Wagner et al., 2018
<https://www.biorxiv.org/content/early/2018/04/09/217737>), but any
UMI-based expression matrix should work equally.

1.  RCA_reference_construction.R # Build references for cell type
    identification
2.  RCA_cellType_identification.R # Cell type identification using
    reference component analysis

### Differnetial expression analysis

To prepare the data from cell type analysis for differential expression
analysis, we run pre-DEG_analysis.R. This script will ensure that the
input to Monocle_v3_RAdata.r is is the correct format. It will also
compute some basic statistics, such as the number of cells per sample
(cell type, time point, etc).

Monocle_v3_RAdata.r is run to calculate differentially expressed genes
(DEGs) between allergen challenged and non-challenged samples. Some
modifications to the code are needed to compute DEGs between any other
pair of groups.

The fold changes (FCs) for the DEGs, we calculate using
sc_FC_zero-infl-neg-bionmial.R. This will compute FCs between allergen
challenged vs non-challenged, meaning that a positive FC indicates a
higher mean expression in the allergen challenged compared to the
non-challenged group. Some modifications to the code are needed to
compute the FCs between any other pair of groups.

We plot the distrinbution of DEGs over different time points and cell
types (Fig 1), as well as the number of cell types over different time
points and treatment groups (Fig 2), by plot_DEGs.R.

1.  pre-DEG_analysis.R # statistics and preparation of input files for
    DEG analysis  
2.  Monocle_v3_RAdata.r # Calculate DEGs
3.  sc_FC_zero-infl-neg-bionmial.R # FC calculations
4.  plot_DEGs.R # Plot the distribution of DEGs

<img src="https://user-images.githubusercontent.com/51739216/144604687-54ee4a7c-b661-4ec7-927b-49d54efe2883.png" width="500" />  
Fig 1. log(number of DEGs) identified between allergen-stimulated and
diluent-stimulated cells, for each cell type at the different time
points.

<img src="https://user-images.githubusercontent.com/51739216/144604759-9677b860-46a9-4296-876b-eb5de3b69b1a.png" width="500" />  
Fig 2. Cell type proportions in different groups of allergen-stimulated
(A) and diluent-stimulated (D), as well as unstimulated control (C),
samples from healthy individuals at the different time points.

## Upstream Regulatory (UR) gene prediction in Ingenuity Pathway Analysis (IPA)

IPA is a commercial software, but you can request a free trial here
(<https://digitalinsights.qiagen.com/products-overview/discovery-insights-portfolio/analysis-and-visualization/qiagen-ipa/>).

To reproduce the analyses for prediction of URs from our project, follow
the following pipeline.

1.  Creat a project (eg. SAR) in the Project Manager to upload the DEGs
    (differentially expressed genes)
    <img src="https://user-images.githubusercontent.com/51739216/155988719-ef25e83b-24ca-4e51-bd9f-d9a4c15de186.png" width="300" />

For each list of DEGs, perform step 2 - 8, for UR prediction in IPA. In
those cases where \>5000 significant DEGs were identified, the DEGs need
to be prioritized, due to limitations in IPA. We included the top 5000
DEGs (based on lowest q-value) into the IPA analysis.

2.  Upload the DEGs into the project “Dataset Files”, including their
    corresponding LogFCs and q-values. Based on this data, choose the
    corresponding ID, “Human gene symbol”, and the observation names,
    “Expr Log Ratio” (LogFC) and “Expr False Discovery Rate” (q-val).
    Keep both LogFC and q-val as the same group, “observation 1”. ”Save”
    and name the dataset.
    <img src="https://user-images.githubusercontent.com/51739216/155989006-ebb0b1b1-e2bb-4677-93b4-63a636739b6e.png" width="450" />

3.  In the lower right corner, click “Analyze/Filter Dataset” and then
    “Core Analysis” to perform IPA analysis of the data.
    <img src="https://user-images.githubusercontent.com/51739216/155743530-21a556a6-26d2-4e3c-95eb-de632dccec07.png" width="450" />

4.  Click ”next”, to get to the settings.
    <img src="https://user-images.githubusercontent.com/51739216/155743696-ac146033-b553-44b1-a020-9ca9e8bae46d.png" width="350" />

5.  In the settings. based on this dataset, define “General settings -
    Species” = Human, “Node Types” = All, “Data Sources” = All,
    “Tissues&Cell Lines” = All, and “Mutation” = All.

6.  Run the analyses by “Run Analysis”.
    <img src="https://user-images.githubusercontent.com/51739216/156144451-22be3f42-08d2-4698-abef-22737098cfe8.png" width="450" />

7.  All the performed analyses can be found in the “SAR” project under
    “Analyses”. Choose one of the analyses to check the results.
    <img src="https://user-images.githubusercontent.com/51739216/156144896-1f4c9511-02df-474c-a4f1-ba39ec1f67ca.png" width="300" />

8.  In the top tab tools, go to “Upstream Analysis” to show the UR
    prediction results. The results can be downloaded by clicking
    <img src="https://user-images.githubusercontent.com/51739216/155744357-86857f73-8111-497d-a716-983ef5abe525.png" width="30" />
    <img src="https://user-images.githubusercontent.com/51739216/155744562-c6748942-bb6f-4fb7-b602-1b0be72e8aaf.png" width="450" />

## Construction of multicellular network models (MNMs) and ranking of URs

To run the whole pipeline, on the example data-set, after IPA UR
prediction, run post_IPA.sh

> ./post_IPA.sh

### MNM construction

The MNMs are constructed for each time point by running
MNM_construction.R in R 4.0. Input to this script are the UR predictions
from IPA and the list of DEGs. For the script to run smoothly, ensure
that the input directory comtaining the UR predictions from IPA includes
one subdirectory per time point, which in turn should include only, but
all, the files to construct the MNMs. Additionally, ensure that the
directory containing the DEGs includes only, but all, the files for MNM
construction. See ‘Rscript MNM_construction.R –help’ for more details.

> Rscript MNM_construction.R ../example_data/IPA_UR-prediction
> ../example_data/DEGs ../output

<img src="https://user-images.githubusercontent.com/51739216/156147705-376f2547-e328-4a4c-9505-9c416da3aa36.png" width="500" />
Fig 3. eg. MNM at 0h, illustrated using cytoscape
(<https://cytoscape.org/>).

### Ranking of URs

The URs from the IPA predictions are ranked based on the number of cell
types and time points in which they were predicted, by running
UR_ranking.R in R 4.0. Input to this script are the UR predictions from
IPA. The structure of the data should be the same as for MNM
construction, with one subdirectory per time point. See ‘Rscript
UR_ranking.R –help’ for more details.

> Rscript UR_ranking.R ../example_data/IPA_UR-prediction ../output

<img src="https://user-images.githubusercontent.com/51739216/156146083-8c047e4f-f135-4b48-b9c9-d8ee52d10fcd.png" width="500" />

Fig 4. Heatmap with rank-ordered URs in rows and all the cell types over
all time points in columns.
