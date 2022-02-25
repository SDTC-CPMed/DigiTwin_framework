module add R/4.0.4

# Create MCDM models
Rscript MNM_construction.R ../example_data/IPA_UR-prediction ../example_data/DEGs ../output

# rank URs
Rscript UR_ranking.R ../example_data/IPA_UR-prediction ../output
