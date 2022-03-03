dir.create("data/Multicellular_Network_Models")
dir.create("data/UR-rank")


### UR prioritization
source('MNM_construction.R')
MNM_construction('data/IPA_UR-prediction', 
                 'data/DEGS_with_Monocle/Monocle_out_withFCs', 
                 'data/Multicellular_Network_Models')
# MNM_construction('data/IPA_UR-prediction', 
#                  'data/DEGS_with_Monocle/Monocle_out_withFCs', 
#                  'data/Multicellular_Network_Models',
#                  time_points = '3D')

### UR prioritization
source('UR_ranking.R')
UR_ranking('data/IPA_UR-prediction', 'data/UR-rank')


# 
# ## Test function
# ### UR prioritization
# source('scripts/MNM_construction.R')
# MNM_construction('IPA_output/IPA_output', 
#                  '/data/sharedData/SAR_allergen_challenge_timeseries/SAR_allergic_vs_healthy_study/MONOCLE/data/DEGS_with_Monocle/Monocle_out/withFCs_timepointX_Sick_vs_Healthy', 
#                  'example_data/test')
# 
# ### UR prioritization
# source('scripts/UR_ranking.R')
# UR_ranking('IPA_output/IPA_output', 'example_data/test')

