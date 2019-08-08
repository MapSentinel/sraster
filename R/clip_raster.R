
#====================
#geotranformation
#====================

geotransform <- function(geo_info,shape){
  m_geo_scale = matrix(c(geo_info[6],0,0, -geo_info[7]),ncol = 2, nrow = 2)
  lefuppercoord_x = geo_info[4]
  lefuppercoord_y = geo_info[5] + geo_info[6]*geo_info[1]
  m_geo_origin = matrix(c(lefuppercoord_x, lefuppercoord_y), nrow = 2, ncol = 1)

  coord_shape_leftupper = matrix(c(st_bbox(shape)[1],st_bbox(shape)[4]), nrow = 2, ncol = 1)
  coord_shape_rightlower = matrix(c(st_bbox(shape)[3],st_bbox(shape)[2]), nrow = 2, ncol = 1)

  m_geo_origin_trans = solve(m_geo_scale) %*% m_geo_origin
  pixel_left_upper = as.integer((solve(m_geo_scale)%*%coord_shape_leftupper) - m_geo_origin_trans)
  pixel_right_lower = as.integer((solve(m_geo_scale)%*%coord_shape_rightlower) - m_geo_origin_trans)

  xcount = abs(pixel_left_upper[1] - pixel_right_lower[1]) + 1
  ycount = abs(pixel_left_upper[2] - pixel_right_lower[2]) + 1

  #coord for clip raster
  coord_pixel_left_upper = m_geo_scale %*% pixel_left_upper + m_geo_origin
  x_max = coord_pixel_left_upper[1] + (xcount)*geo_info[6]
  y_min = coord_pixel_left_upper[2] - (ycount)*geo_info[7]
  bbox_sraster = c(coord_pixel_left_upper[1], y_min, x_max, coord_pixel_left_upper[2])
  names(bbox_sraster) = c("xmin","ymin","xmax","ymax")

  result = list()
  result[["offs"]] = c(pixel_left_upper[2], pixel_left_upper[1])
  result[["xcount"]] = xcount
  result[["ycount"]] = ycount
  result[["bbox"]] = bbox_sraster

  return(result)
}


#===============================
#Importing image
#===============================
#For this example we call only one image
library(rgdal)

read_raster <-function(files_raster, xyoffs, reg, r_shape_array, shape){
  list_bands = list()
  n_img = length(files_raster)
  pb = txtProgressBar(min = 1, max = n_img,style = 3)
  for(j in 1:length(files_raster)){
    geo_info = rgdal::GDALinfo(files_raster[j],returnStats = FALSE)
    datasource = rgdal::GDAL.open(files_raster[j])
    #This part needs to be improved

    for(i in 1:geo_info[3]){
      band = paste0('IMG_',j,"_B_",i)
      clip_extent = rgdal::getRasterData(datasource, offset = xyoffs, region.dim = reg, band=i)
      if(i == 2){
        red_band = clip_extent * r_shape_array
      }
      else if(i == 3){
        near_band = clip_extent * r_shape_array
      }
      list_bands[[band]] = clip_extent * r_shape_array
    }
    #NDVI
    ndvi = (near_band - red_band)/(near_band + red_band)
    name_ndvi = paste0("NDVI_",j)
    list_bands[[name_ndvi]] = ndvi
    rgdal::closeDataset(datasource)
    red_band = 0
    near_band = 0
    Sys.sleep(0.1)
    setTxtProgressBar(pb,j)
  }
  multiarray = array(unlist(list_bands),dim = c(reg[2], reg[1] , length(list_bands)))
  result = list()
  result[["bands"]] <- names(list_bands)
  result[["array"]] <- multiarray
  return(result)
}


#crating sraster function
as_sraster <- function(shape,files_raster){
  cat(paste0("\tShape ",shape$OBJECTID),"\n")
  geo_info = rgdal::GDALinfo(files_raster[1],returnStats = FALSE)
  coord_pixel = geotransform(geo_info, shape)
  xyoffs = coord_pixel$offs
  nrows = coord_pixel$ycount
  ncols = coord_pixel$xcount
  reg = c(nrows,ncols)

  #Rasterize
  r <- stars::st_as_stars(st_bbox(coord_pixel$bbox) , nx = ncols, ny = nrows )
  shape$field = 1
  r_shape = stars::st_rasterize(shape[,"field"], r,options = "ALL_TOUCHED=TRUE")
  r_shape$field[r_shape$field == 0] = NA
  r_shape_array = r_shape$field

  result <- read_raster(files_raster, xyoffs, reg, r_shape_array,shape)

  #Coordinates

  x_min = coord_pixel$bbox[1] + geo_info[6]/2
  y_min = coord_pixel$bbox[2] + geo_info[7]/2
  x_max = coord_pixel$bbox[3] - geo_info[6]/2
  y_max = coord_pixel$bbox[4] - geo_info[7]/2

  x_line = seq(x_min,x_max,geo_info[6])
  y_line = seq(y_max,y_min,-geo_info[7])
  coord_x_matrix = array(rep(x_line, ncols),dim = c(ncols, nrows))
  coord_y_matrix = t(array(rep(y_line, nrows),dim = c(nrows, ncols)))

  coordinates = data.frame(x=c(coord_x_matrix), y=c(coord_y_matrix))
  result$coordinates = coordinates

  class(result) <- "sraster"
  return(result)
}


print.sraster <-
  function(x,...){
    cat("class:    sraster", "\n")
    cat("space dimension: nrows: ", dim(x$array)[1], " nccols: ", dim(x$array)[2],"\n")
    cat("Number of layer: ", dim(x$array)[3],"\n")
  }

#===============================
#Clustering
#===============================
library(reshape2)

convert_2d <- function(x){
  y = reshape2::melt(x)
  return(y[,3])
}

func_nan <- function(x){all(is.na(x))}

library(fpc)

kmeans_sraster <-
  function(x,...)
  {
    array_sr = x$array

    data_2d = apply(array_sr,3,convert_2d)

    index_nonan = apply(data_2d,1,func_nan)
    data_2d_nonan = data_2d[!index_nonan,]
    #I need to improve the following
    data_2d_nonan[is.na(data_2d_nonan)] <- 999999

    #defining number of clusters
    a = c(1:5)
    for(l in 2:6){
      km <- kmeans(data_2d_nonan,l)
      a[l-1] = round(fpc::calinhara(data_2d_nonan,km$cluster),digits=2)
    }
    centers = order(a)[1]
    cluster_data_2d = kmeans(data_2d_nonan,centers)

    vector_cluster = c(array(NA, dim = dim(array_sr)[1:2]))

    vector_cluster[!index_nonan]<-cluster_data_2d$cluster

    matrix_cluster = array(vector_cluster, dim = dim(array_sr)[1:2])
    return(matrix_cluster)
  }



#===============================
#Clip sraster
#===============================

clip_sraster = function(x, mask, type){
  if(type == 'Majority rule'){
    largest_cluster = order(table(mask),decreasing = TRUE)[1]
    mask[mask != largest_cluster] <- NA
    mask[mask == largest_cluster] <- 1
  }
  else if(type == 'rule_ndvi'){
    bands_ndvi = grep("NDVI",x$bands)
    cluster_n = as.numeric(names(table(mask)))
    median_ndvi = c(1:length(cluster_n)) 
    median_ndvi[]<-0
    for(i in cluster_n){
      list_ndvi = list()
      for(j in bands_ndvi){
          bndvi = x$array[,,j]
          list_ndvi[[j]] = c(bndvi[mask==i])
      }
      median_ndvi[i] = median(unlist(list_ndvi),na.rm = TRUE)
    }
    lowest_ndvi = order(median_ndvi)[1]
    mask[mask != lowest_ndvi] <- NA
    mask[mask == lowest_ndvi] <- 1
  }

  func_mask <- function(y, mask) y * mask

  array_clip_list = lapply(asplit(x$array,3),func_mask,mask)
  array_clip = array(unlist(array_clip_list),dim = dim(x$array))
  x$array <- array_clip
  return(x)
}

#===============================
#as data frame sf
#===============================

as.data.frame.sraster <-
  function(x,...)
  {
    array_x = x$array
    data_2d = apply(array_x,3,convert_2d)
    index_nonan = apply(data_2d,1,func_nan)

    data_2d_nonan = data_2d[!index_nonan,]
    coordxy = x$coordinates[!index_nonan,]
    df_data_2d_nonan = data.frame(coordxy, data_2d_nonan)
    colnames(df_data_2d_nonan) <- c("x","y",x$bands)

    df_spatial = st_as_sf(df_data_2d_nonan, coords = c("x","y"))
    return(df_spatial)
  }


#===============================
#Example
#===============================

library(rgdal)
library(sf)
library(stars)
library(raster)
#setwd("C://IPSTERS")
#file_shape = 'C:\\IPSTERS\\COS2015\\COS2015_v2_08_02_2019_clip.shp'
file_shape = '/home/willimarti2008/Documents/DGT/COS2015/COS2015_v2_08_02_2019_clip.shp'
legend = st_read(file_shape)
legend_water = legend[which(legend$Legend == 'water'),]

#===============================
#random selection
#===============================
#Goal : stratified random selection of traing samples at level of polygon, (queriying only one class)

n_samples = 20
set.seed(123)   #setting same random selection for testing
index = sample(1:dim(legend_water)[1], n_samples,replace = FALSE)
query_water = legend_water[index,]
#===============================
#Calling imagery
#===============================

#images_folder = 'C:\\IPSTERS\\IMAGES'
images_folder = '/home/willimarti2008/Documents/DGT/images'
join_path = function(x,path_folder){
  a = strsplit(x,'[.]')
  format_file = a[[1]][length(a[[1]])]
  if(format_file == 'tif'){
    return(paste0(path_folder,'/',x)) 
  }
}

list_images1 = list.files(images_folder)
list_images2 = c("S2A_L2A_20171002-113001_T29SND.tif",
                 "S2A_L2A_20171121-112837_T29SND.tif",
                 "S2A_L2A_20171221-112810_T29SND.tif",
                 "S2A_L2A_20180321-112321_T29SND.tif",
                 "S2A_L2A_20180619-112602_T29SND.tif",
                 "S2A_L2A_20180729-112845_T29SND.tif",
                 "S2A_L2A_20180818-112627_T29SND.tif",
                 "S2A_L2A_20180927-112959_T29SND.tif",
                 "S2A_L2A_20181007-112305_T29SND.tif")



paths_images = unlist(lapply(list_images1,join_path,images_folder))

list_shapes = split(query_water,query_water$OBJECTID)

shape = query_water[1,]

workflow <- function(shape,paths_images){
          result = as_sraster(shape,paths_images)
          #plot(raster((result$array)[,,1]))
          cluster_matrix = kmeans_sraster(result)
          #plot(raster(cluster_matrix))
          result_clip = clip_sraster(result, mask = cluster_matrix, type ="Majority rule")
          #plot(raster((result_clip$array)[,,1]))
          return(as.data.frame(result_clip))
          }

result_list = lapply(list_shapes,workflow,paths_images)

result_df = do.call("rbind",result_list)
plot(st_geometry(result_df))


st_write(result_df, "output.csv", layer_options = "GEOMETRY=AS_XY")



