
#====================
#geotranformation
#====================

geotransform <- function(geo_info,shape){
  m_geo_scale = matrix(c(geo_info[6],0,0, -geo_info[7]),ncol = 2, nrow = 2)
  lefuppercoord_x = geo_info[4]
  lefuppercoord_y = geo_info[5] + geo_info[6]*geo_info[1]
  m_geo_origin = matrix(c(lefuppercoord_x, lefuppercoord_y), nrow = 2, ncol = 1)

  coord_shape_leftupper = matrix(c(extent(shape)[1],extent(shape)[4]), nrow = 2, ncol = 1)
  coord_shape_rightlower = matrix(c(extent(shape)[2],extent(shape)[3]), nrow = 2, ncol = 1)

  m_geo_origin_trans = solve(m_geo_scale) %*% m_geo_origin
  coord_pixel_left_upper = as.integer((solve(m_geo_scale)%*%coord_shape_leftupper) - m_geo_origin_trans)
  coord_pixel_right_lower = as.integer((solve(m_geo_scale)%*%coord_shape_rightlower) - m_geo_origin_trans)
  result = list()
  result[["offs"]] = c(coord_pixel_left_upper[2], coord_pixel_left_upper[1])
  result[["xcount"]] = abs(coord_pixel_left_upper[1] - coord_pixel_right_lower[1])
  result[["ycount"]] = abs(coord_pixel_left_upper[2] - coord_pixel_right_lower[2])
  return(result)
}


#===============================
#Importing image
#===============================
#For this example we call only one image
library(rgdal)
library(raster)

read_raster <-function(files_raster, xyoffs, reg){
  list_bands = list()

  for(j in 1:length(files_raster)){
    geo_info = rgdal::GDALinfo(files_raster[j],returnStats = FALSE)
    datasource = rgdal::GDAL.open(files_raster[j])
    for(i in 1:geo_info[3]){
      band = paste0('IMG_',j,"_B_",i)
      list_bands[[band]] = rgdal::getRasterData(datasource, offset = xyoffs, region.dim = reg, band=i)
    }
    rgdal::closeDataset(datasource)
  }
  multiarray = array(unlist(list_bands),dim = c(reg[2], reg[1] , length(list_bands)))
  result = list()
  result[["bands"]] <- names(list_bands)
  result[["array"]] <- multiarray
  return(result)
}



#example
crop_raster <- function(shape,files_raster){
  geo_info = rgdal::GDALinfo(files_raster[1],returnStats = FALSE)
  coord_pixel = geotransform(geo_info, shape)
  xyoffs = coord_pixel$offs
  nrows = coord_pixel$ycount
  ncols = coord_pixel$xcount
  reg = c(nrows,ncols)
  result <- read_raster(files_raster, xyoffs, reg)
  return(result)
}


#===============================
#Example
#===============================


library(rgdal)
#setwd("C://IPSTERS")
file_shape = 'C:\\IPSTERS\\COS2015\\COS2015_v2_08_02_2019_clip.shp'
#legend = st_read(file_shape)


#===============================
#random selection
#===============================
#Goal : stratified random selection of traing samples at level of polygon, (queriying only one class)


legend = rgdal::readOGR(file_shape)
n_samples = 30
list_classes = split(legend,legend$Legend)
set.seed(123)   #setting same random selection for testing
index = sample(1:length(list_water$water), n_samples,replace = FALSE)
query_water = list_classes$water[index,]

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


shape = split(query_water,query_water$OBJECTID)

result = lapply(shape,crop_raster,files_raster[c(1:12)])



#===============================
#Rasterization
#===============================

#rasterize
#r <- raster(ncols = reg[1], nrows = reg[2])
#extent(r) <- extent(shape)
#r <- rasterize(shape, r)


