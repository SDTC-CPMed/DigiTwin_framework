dir.create("plot")
dir.create("data")
dir.create("data/DEGS_with_Monocle")
dir.create("data/DEGS_with_Monocle/Matrix_in")
dir.create("data/DEGS_with_Monocle/Monocle_out")
dir.create("data/DEGS_with_Monocle/Monocle_out_withFCs")
dir.create("data/knn_smoothing")
dir.create("data/RCA_out")
dir.create("data/RCA_out/full_matrix")
dir.create("data/RCA_references")


source('sc_data_quality_sorting.R')
source('RCA_reference_construction.R')
source('RCA_cellType_identification.R')
source('pre-DEG_analysis.R')
source('Monocle_v3_RAdata.r')
source('sc_FC_zero-infl-neg-bionmial.R')
source('plot_DEGs.R')



###
# Read data from GSE180697 and save it to data folder
#filenames = c('data/GSE180697_SAR_patients_expression_matrix.csv.gz',
#              'data/GSE180697_Healthy_controls_expression_matrix.csv.gz')

filenames = c('/data/sharedData/SAR_allergen_challenge_timeseries/SAR_allergic_patient_study/DGE_data/HA_min200genesPerCell_UMI_expression_matrix.txt.gz',
              '/data/sharedData/SAR_allergen_challenge_timeseries/SAR_healty_control_study/DGE_data/HC_min200genesPerCell_UMI_expression_matrix.txt.gz')

HA = sc_data_quality_sorting(filenames[1])
write.csv(HA, 'data/HA_min200genesPerCell_sorted_expression_matrix.csv', quote = F)
HC = sc_data_quality_sorting(filenames[2])
write.csv(HC, 'data/HC_min200genesPerCell_sorted_expression_matrix.csv', quote = F)



### Run knn smoothing
#system(chmod u+x run_knn_smoothing.sh)

system('./run_knn_smoothing.sh ')


dir.create("plot")
dir.create("data")
dir.create("data/DEGS_with_Monocle")
dir.create("data/DEGS_with_Monocle/Matrix_in")
dir.create("data/DEGS_with_Monocle/Monocle_out")
dir.create("data/DEGS_with_Monocle/Monocle_out_withFCs")
dir.create("data/knn_smoothing")
dir.create("data/RCA_out")
dir.create("data/RCA_out/full_matrix")
dir.create("data/RCA_references")

###
# Read data from GSE180697 and save it to data folder




### Run cell typing
RCA_cellType_identification('HA')
#RCA_cellType_identification('HC')


### Run pre-DEG analysis
pre_DEG_analysis('HA')
#pre_DEG_analysis('HC')

### Run Monocle_v3_RAdata
Monocle_v3_RAdata('HA')
#Monocle_v3_RAdata('HC')

### Run sc_FC_zero_infl_neg-binomial
sc_FC_zero_infl_neg_binomial('HA')
#sc_FC_zero_infl_neg_binomial('HC')

### Run plot_DEGs.R
plot_custom()
