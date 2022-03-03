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
# Code by https://github.com/yanailab/knn-smoothing
system('./run_knn_smoothing.sh ')



### Run cell typing
source('RCA_reference_construction.R')
source('RCA_cellType_identification.R')
RCA_cellType_identification('HA')
RCA_cellType_identification('HC')

### Run pre-DEG analysis
source('pre-DEG_analysis.R')
pre_DEG_analysis('HA')
pre_DEG_analysis('HC')

### Run Monocle_v3_RAdata
source('Monocle_v3_RAdata.r')

#Type of analysis could be 'HA_vs_HC' for DEGs between healthy and sick,
#and 'HA' or 'HC' for dilutant vs allergen challenged
Type_of_Analyses = 'HA_vs_HC'
Monocle_v3_RAdata(Type_of_Analyses)



### Run sc_FC_zero_infl_neg-binomial
source('sc_FC_zero-infl-neg-bionmial.R')
sc_FC_zero_infl_neg_binomial(Type_of_Analyses)



### Run plot_DEGs.R
source('plot_DEGs.R')
plots = plot_custom()
pdf('plot/lognDEGs_over_time.pdf')
plots[[1]]
dev.off()
pdf('plot/HA_celltype_ratios_over_groups.pdf')
plots[[2]]
dev.off()
pdf('plot/HC_celltype_ratios_over_groups.pdf')
plots[[3]]
dev.off()


