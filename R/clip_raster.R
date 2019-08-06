
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

    for(i in 1:geo_info[3]){
      band = paste0('IMG_',j,"_B_",i)
      clip_extent = rgdal::getRasterData(datasource, offset = xyoffs, region.dim = reg, band=i)
      list_bands[[band]] = clip_extent * r_shape_array
    }
    rgdal::closeDataset(datasource)
    Sys.sleep(0.1)
    setTxtProgressBar(pb,j)
  }
  multiarray = array(unlist(list_bands),dim = c(reg[2], reg[1] , length(list_bands)))
  result = list()
  result[["bands"]] <- names(list_bands)
  result[["array"]] <- multiarray
  return(result)
}


#example
crop_raster <- function(shape,files_raster){
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
  return(result)
}


#===============================
#Example
#===============================


library(rgdal)
library(sf)
library(stars)
library(raster)
#setwd("C://IPSTERS")
file_shape = 'C:\\IPSTERS\\COS2015\\COS2015_v2_08_02_2019_clip.shp'
legend = st_read(file_shape)
legend_water = legend[which(legend$Legend == 'water'),]

#===============================
#random selection
#===============================
#Goal : stratified random selection of traing samples at level of polygon, (queriying only one class)

n_samples = 5
set.seed(123)   #setting same random selection for testing
index = sample(1:dim(legend_water)[1], n_samples,replace = FALSE)
query_water = legend_water[index,]
#===============================
#Calling imagery
#===============================

images_folder = 'C:\\IPSTERS\\IMAGES'
join_path = function(x,path_folder){
  a = strsplit(x,'[.]')
  format_file = a[[1]][length(a[[1]])]
  if(format_file == 'tif'){
    return(paste0(path_folder,'\\',x))
  }
}

files_raster = unlist(lapply(list.files(images_folder),join_path,images_folder))

list_shapes =split(query_water,query_water$OBJECTID)

result = lapply(list_shapes,crop_raster,files_raster[c(1:12)])

shape = query_water[query_water$OBJECTID == 91801,]


#===============================
#Rasterization
#===============================



