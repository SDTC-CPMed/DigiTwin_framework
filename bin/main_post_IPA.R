dir.create("data/Multicellular_Network_Models")
dir.create("data/UR-rank")


### MNM construction
source('MNM_construction.R')
MNM_construction('data/IPA_UR-prediction', 
                 'data/DEGS_with_Monocle/Monocle_out_withFCs', 
                 'data/Multicellular_Network_Models')
### MNM construction, specific time point
# MNM_construction('data/IPA_UR-prediction', 
#                  'data/DEGS_with_Monocle/Monocle_out_withFCs', 
#                  'data/Multicellular_Network_Models',
#                  time_points = '3D')


### UR prioritization
source('UR_ranking.R')
UR_ranking('data/IPA_UR-prediction', 'data/UR-rank')


