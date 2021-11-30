#' Generate Reference Component featrues
#' 
#' Compute Reference Component features for clustering analysis
#' 
#' @param obj_in data object.
#' @param method choose from three RCA options: "GlobalPanel"(default), "ColonEpitheliumPanel", or "SelfProjection". 
#' @param power power to raise up to for the RCA features before clustering. Default is 4. 
#' @return data object
#' @export
#' @examples
#' 
#' data_obj = featureConstruct(data_obj);
#' 
#' 
featureConstruct <- function(obj_in, method="GlobalPanel", sysdata=NULL, power=4)
{
  require(picante)
#### 1: reading input ####
  if(is.null(sysdata)){
      data(sysdata, envir=environment())
  }
  fpkm_temp = obj_in$fpkm_transformed;
  
  #### 2: choosing method ####
 if (method == "GlobalPanel"){
        data=fpkm_temp;
        data11=list();
        for (i in 1:length(sysdata[[1]])){
        data1 = sysdata[[1]][[i]];
        d1 = as.dist(1-cor(data1));
        t1 = hclust(d1,method="average");
        temp = rownames(data);
        temp1 = gsub("^.*?_","",temp);
        temp2 <- strsplit(temp1,"_ENS");
        temp3 <- paste("",lapply(temp2,"[[",1),sep="");
        temp4 = intersect(temp3,rownames(data1));
        temp5 = temp[match(temp4,temp3)];
        data2 = data1[temp4,];
        data4 = data[temp5,,drop=FALSE];
        data3 = data2;
        data3[data3<=(sysdata$at)[i]] = (sysdata$at)[i];
        data5 = cbind(data4,data3);
        temp = cor.table(data5,cor.method = "pearson");
        data6 = temp$r
        #data6 = cor(data5,method = "pearson");
        data6p = temp$P;
        data6 = as.data.frame(data6);
        data6p = as.data.frame(data6p);
        data7 = data6[(dim(fpkm_temp)[2]+1):dim(data6)[2],1:dim(fpkm_temp)[2]];
        #data7 = data6[(dim(data)[2]+1):dim(data6)[2],1:dim(data)[2]];
        data7p = data6p[(dim(fpkm_temp)[2]+1):dim(data6p)[2],1:dim(fpkm_temp)[2]];
        data8 = data7;
        #data9 = abs(data8)^(sysdata$sp) * sign(data8);
        data9 = abs(data8)^(power) * sign(data8);
        data10 = scale(data9,center=TRUE,scale=TRUE);
        #data11[[i]]=data10;
    }
        #data12 = as.data.frame(t(cbind(as.data.frame(t(data11[[1]])),as.data.frame(t(data11[[2]])))));
        #fpkm_for_clust = data12;
        fpkm_for_clust = data10;
 }
  if (method == "ColonCellLinesGSE10843_avg_multiple"){
	data11=list();	
	data0 = sysdata[[which(names(sysdata)==method)]]
	fpkm_for_clust=list()
	for (z in 1: length(data0)){
		data1 = data0[[z]]
		#if (method == "ColonCellLinesGSE10843_avg"){
		#	data1 = sysdata[[8]]; #ColonCellLinesGSE10843_avg
		#}
		#if (method == "CellRef_Huan_avg_cell"){
		#	data1 = sysdata[[7]]; #}
		d1 = as.dist(1-cor(data1));
		t1 = hclust(d1,method="average");
		temp = rownames(fpkm_temp);
		temp4 = intersect(temp,rownames(data1));
		temp5 = temp[match(temp4,temp)];
		data2 = data1[temp4,];
		data4 = fpkm_temp[temp5,,drop=FALSE];
		data3 = data2;
		data3[data3<=1.5] = 1.5;
		data5 = cbind(data4,data3);
		#data5 = data5[,colSums(data5)!=0]
		data6 = cor(data5,method = "pearson");
		data6 = as.data.frame(data6);
		data7 = data6[(dim(fpkm_temp)[2]+1):dim(data6)[2],1:dim(fpkm_temp)[2]];
		#data7 = data6[(dim(data6)[2]+1):dim(data6)[2],1:dim(data6)[2]];
		data8 = data7;
		data9 = abs(data8)^(power) * sign(data8);
		data10 = scale(data9,center=TRUE,scale=TRUE); 
		fpkm_for_clust[[z]] = data10; 
	
	}
  }
 #if (method == "ColonCellLinesGSE10843_avg"){
#	#data=fpkm_temp;
#    data11=list();
	#data1 = sysdata[[8]]; #ColonCellLinesGSE10843_avg
	#d1 = as.dist(1-cor(data1));
	#t1 = hclust(d1,method="average");
	#temp = rownames(fpkm_temp);
	#temp4 = intersect(temp,rownames(data1));
	#temp5 = temp[match(temp4,temp)];
	#data2 = data1[temp4,];
	#data4 = fpkm_temp[temp5,,drop=FALSE];
	#data3 = data2;
	#data3[data3<=1.5] = 1.5;
	#data5 = cbind(data4,data3);
	#data6 = cor(data5,method = "pearson");
	#data6 = as.data.frame(data6);
	#data7 = data6[(dim(fpkm_temp)[2]+1):dim(data6)[2],1:dim(fpkm_temp)[2]];
	#data8 = data7;
	#data9 = abs(data8)^(power) * sign(data8);
	#data10 = scale(data9,center=TRUE,scale=TRUE); 
	#fpkm_for_clust = data10; 
 #}

 
  #if (method = "IMMGEN_strict" || method == "IMMGEN" || method == "ColonCellLinesGSE10843_avg" ||method == "CellRef_Huan_avg_cell_strict" || method =="BioGPS_MOE430_GSE10246_GeneSymbol_upper"||method == "BioGPS_MOE430_GSE10246"|| method == "CellRef_Huan_avg_cell" || method == "CellRef_Huan" || method == "ColonCellLinesGSE10843"){
  if (!(method %in% c('ColonCellLinesGSE10843_avg_multiple','GlobalPanel','ColonEpitheliumPanel','SelfProjection'))){
      cat("Am I in?\n")
  #data=fpkm_temp;
    data11=list();	
      data1 = sysdata[[which(names(sysdata)==method)]]
      ##if (method == "ColonCellLinesGSE10843_avg"){
      ##	data1 = sysdata[[8]]; #ColonCellLinesGSE10843_avg
      ##}
      ##if (method == "CellRef_Huan_avg_cell"){
      ##	data1 = sysdata[[7]]; #}
      d1 = as.dist(1-cor(data1));
      t1 = hclust(d1,method="average");
      temp = rownames(fpkm_temp);
      temp4 = intersect(temp,rownames(data1));
      if(length(temp4) == 0){
          return(NULL)
      }
      temp5 = temp[match(temp4,temp)];
      data2 = data1[temp4,];
      data4 = fpkm_temp[temp5,,drop=FALSE];
      data3 = data2;
      data3[data3<=1.5] = 1.5;
      data5a = cbind(data4,data3);
      data5 = data5a[,colSums(data5a)!=0];
      temp = cor.table(data5,cor.method = "pearson");
      data6 = temp$r;
      data6p = temp$P;
      data6 = as.data.frame(data6);
      data6p = as.data.frame(data6p);
      data7 = data6[(dim(fpkm_temp)[2]+1):dim(data6)[2],1:dim(fpkm_temp)[2]];
      data7p = data6p[(dim(fpkm_temp)[2]+1):dim(data6p)[2],1:dim(fpkm_temp)[2]];
                                        #data7 = data6[(dim(data6)[2]+1):dim(data6)[2],1:dim(data6)[2]];
      data8 = data7;
      data9 = abs(data8)^(power) * sign(data8);
      data10 = scale(data9,center=TRUE,scale=TRUE); 
      fpkm_for_clust = data10; 
 }
 
 if (method == "ColonEpitheliumPanel"){
   
      fc = apply(sysdata$ColonEpiPanel,1,function(x) x-median(x));
      fs = fc>1.5;
      fs1 = rownames(sysdata$ColonEpiPanel[apply(fs,1,function(x) sum(x))>0,]);
      gl_intersect = intersect(rownames(fpkm_temp),fs1);
      data_combine0 = cbind(fpkm_temp[gl_intersect,],sysdata$ColonEpiPanel[gl_intersect,]);
      data_combine = data_combine0;
      cor_matrix = cor(data_combine,method = "pearson");
      cor_matrix = as.data.frame(cor_matrix);
      individual_cell_type_prediction_matrix = cor_matrix[(dim(fpkm_temp)[2]+1):dim(cor_matrix)[2],1:dim(fpkm_temp)[2]];
      fpkm_for_clust_0 = individual_cell_type_prediction_matrix;
      fpkm_for_clust_power = abs(fpkm_for_clust_0)^power * sign(fpkm_for_clust_0);
      fpkm_for_clust = fpkm_for_clust_power;  
 }


 
 if (method == "SelfProjection"){    

    min_val = min(fpkm_temp);
    fuzzy_transfer <- function(x,min_val){
      upper_thresh = max(log10(50),median(x[x>min_val]));
      y=x;
      y[x==min_val]=0;
      y[x>=upper_thresh]=1;
      y[x<upper_thresh & x>min_val]=(x[x<upper_thresh & x>min_val]-min_val)/upper_thresh;      
      return(y);
    }
    
    fpkm_transfer = apply(fpkm_temp,1,function(x) fuzzy_transfer(x,min_val));
    fpkm_transfer1 = as.data.frame(t(fpkm_transfer));
    ones_vec = apply(fpkm_transfer1,1,function(x) sum(x==1));
    zeros_vec = apply(fpkm_transfer1,1,function(x) sum(x==0));

    fs_vec_fuzzy = ones_vec >= dim(fpkm_temp)[2]*1/15 & 
                   zeros_vec >= dim(fpkm_temp)[2]*1/15 ;    
    fs_vec = fs_vec_fuzzy;
    fpkm_for_clust0 = fpkm_temp[fs_vec,,drop=FALSE];
    
    fpkm_for_clust = as.data.frame(cor(fpkm_for_clust0,method="pearson"));
  }

  #### 3: writing output ####  
  obj_out = append(obj_in,
                list("fpkm_for_clust" = fpkm_for_clust,
					"corrScorP" = data7p,
					"corrScorRHO" = data7)
  );
  return(obj_out)
}
