#' Generate cell clusters
#' 
#' Hierarchical clustering and dynamic tree cutting. 
#' 
#' @param obj_in data object.
#' @param method can only be "hclust" (default). 
#' @param deepSplit_wgcna integer value indicating how deep the dendrogram should be cut. Values can be 0, 1 (default), 2, 3 and 4. 
#' @param min_group_Size_wgcna integer value indicating the minimum size of the resulting clusters. Default is 5.
#' @return data object.
#' @export
#' @examples
#' 
#' data_obj = cellClust(data_obj);
#' 
cellClust <- function(obj_in,method="hclust",deepSplit_wgcna=1,min_group_Size_wgcna=5)
{
  #### 1: reading input ####  
  fpkm_temp = obj_in$fpkm_for_clust;
  
  if (!require(flashClust)) install.packages("flashClust",repos = "http://cran.us.r-project.org") 
  require(flashClust)
  if (!require(WGCNA)){
    source("http://bioconductor.org/biocLite.R");
    biocLite(c("impute", "GO.db", "preprocessCore")); 
    install.packages("WGCNA");
  }
  require(WGCNA)  
  
  #### 2: choosing method ####
    if (method == "hclust"){
		if (!is.list(fpkm_temp)){
			d = as.dist(1-cor(fpkm_temp,method="pearson"));
			cellTree = flashClust(d,method = "average");     
			dynamicGroups = cutreeDynamic(dendro = cellTree,distM = as.matrix(d),deepSplit = deepSplit_wgcna,
										  pamStage = FALSE,minClusterSize= min_group_Size_wgcna);
			dynamicColors = labels2colors(dynamicGroups);
			if (!is.null(names(fpkm_temp))){
					group_labels = matrix(dynamicGroups,nrow = length(names(fpkm_temp)),ncol = 1,list(names(fpkm_temp),c("groupLabel")),byrow=FALSE)
				}else{
					group_labels = matrix(dynamicGroups,nrow = length(colnames(fpkm_temp)),ncol = 1,list(colnames(fpkm_temp),c("groupLabel")),byrow=FALSE)
				}
			group_labels = as.data.frame(group_labels)
			group_labels_color = cbind(group_labels,dynamicColors);
		}else{
			group_labels = list()
			group_labels_color = list()
			cellTree = list()
			dynamicColors=list()
			d = list()
			for (z in 1:length(fpkm_temp)){
				d[[z]] = as.dist(1-cor(fpkm_temp[[z]],method="pearson"));
				cellTree[[z]] = flashClust(d[[z]],method = "average");     
				dynamicGroups = cutreeDynamic(dendro = cellTree[[z]],distM = as.matrix(d[[z]]),deepSplit = deepSplit_wgcna,
											  pamStage = FALSE,minClusterSize= min_group_Size_wgcna);
				dynamicColors[[z]] = labels2colors(dynamicGroups);
				if (!is.null(names(fpkm_temp[[z]]))){
						group_labels[[z]] = matrix(dynamicGroups,nrow = length(names(fpkm_temp[[z]])),ncol = 1,list(names(fpkm_temp[[z]]),c("groupLabel")),byrow=FALSE)
					}else{
						group_labels[[z]] = matrix(dynamicGroups,nrow = length(colnames(fpkm_temp[[z]])),ncol = 1,list(colnames(fpkm_temp[[z]]),c("groupLabel")),byrow=FALSE)
					}
				group_labels[[z]] = as.data.frame(group_labels[[z]])
				group_labels_color[[z]] = cbind(group_labels[[z]],dynamicColors[[z]]);
			}
		}
    }
  #### 3: writing output ####    
  obj_out = append(obj_in,
              list("d" = d,
              "cellTree" = cellTree,
              "dynamicColors" = dynamicColors,
              "group_labels_color" = group_labels_color
              )
  )
  
  return(obj_out)
  
}