


file_directory = 'C:\\IPSTERS\\sraster\\ins4\\sampling_01'
list_files = list.files(file_directory)
file_directory_output = 'C:\\IPSTERS\\sraster\\ins4\\sampling_01'


for(file_x in list_files){
  
  x = read.csv(paste0(file_directory,"\\",file_x),sep = ",", header = TRUE)
  names_x = colnames(x)
  ind_NDVI = grep("NDVI",names_x)
  ind_NDBI = grep("NDBI",names_x)
  ind_NDMIR = grep("NDMIR",names_x)
  
  indeces = list(x[,ind_NDVI],x[,ind_NDBI],x[,ind_NDMIR])
  indeces_names = c("NDVI","NDBI","NDMIR")
  
  statistic_mean = function(y) mean(y, na.rm=TRUE)
  statistic_min = function(y) min(y, na.rm=TRUE)
  statistic_max = function(y) max(y, na.rm=TRUE)
  statistic_var = function(y) var(y, na.rm=TRUE)
  statistic_q10 = function(y) quantile(y, probs = c(0.1), na.rm=TRUE)
  statistic_q25 = function(y) quantile(y, probs = c(0.25), na.rm=TRUE)
  statistic_q50 = function(y) quantile(y, probs = c(0.50), na.rm=TRUE)
  statistic_q75 = function(y) quantile(y, probs = c(0.75), na.rm=TRUE)
  statistic_q90 = function(y) quantile(y, probs = c(0.90), na.rm=TRUE)

  statistic = list(statistic_mean, statistic_min, statistic_max, statistic_var, statistic_q10, statistic_q25, statistic_q50, statistic_q75, statistic_q90)
  statistic_names = c("mean","min","max","var","q10","q25","q50","q75","q90")
  output = NULL
  for(j in 1:length(statistic)){
    statistic_w = lapply(indeces, function(w) apply(w,1,statistic[[j]]))
    df_stat = do.call("cbind",statistic_w)
    colnames(df_stat)<- paste0(indeces_names, statistic_names[j])
    output = cbind(output,df_stat)
  }
  output_df = as.data.frame(output)
  
  #merging original dataframe with the new one
  
  x_output = cbind(x,output_df)
  
  #saving file

  #write.csv(x_output,paste0(file_directory_output,"\\",file_x))
  #x = NULL
  #x_output = NULL
  cat('done ',file_x)
}

  
  