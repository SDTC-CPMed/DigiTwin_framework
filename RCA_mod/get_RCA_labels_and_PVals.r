get_RCA_labels_and_PVals = function(data_obj) {

 	 cells = rownames(data_obj$fpkm_for_clust)
	 label = cbind(data_obj$group_labels_color,data.frame(hg = colnames(data_obj$fpkm_for_clust)))
 
 	 for (z in 1:ncol(data_obj$fpkm_for_clust)){
 	   label$label_max_RHO[z] = cells[data_obj$fpkm_for_clust[,z]==max(data_obj$fpkm_for_clust[,z])]
  	   label$Pval_max_RHO[z] = data_obj$corrScorP[data_obj$corrScorRHO[,z]==max(data_obj$corrScorRHO[,z]),z]
 	 }
	 # get RCA label:
	for (z in 1:max(label$groupLabel)){
		label$label_max_RHO[label$groupLabel==z]
		label$RCA[label$groupLabel==z] = names(which.max(table(label$label_max_RHO[label$groupLabel==z])))
	}
	# get p-values:
	for (z in 1:ncol(data_obj$fpkm_for_clust)){
		label$RCA_Pval[z] = data_obj$corrScorP[label$RCA[z],z]
	}
	label$RCA_Pval_adj = p.adjust(label$RCA_Pval,method='BH')
	label
}
