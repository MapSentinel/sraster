#===============================
#Importing shapefile
#===============================

library(sf)
setwd("/home/willimarti2008/Documents/DGT/sraster")
file_shape = 'COS2015/COS2015_v2_08_02_2019_clip.shp'
#legend = st_read(file_shape)

legend = rgdal::readOGR(file_shape)

#===============================
#random selection 
#===============================
#Goal : stratified random selection of traing samples at level of polygon

set.seed(123)   #setting same random selection for testing
list_df_classes  = split(legend,legend$Legend)

strat_random_selection = function(x, n_samples){
  n = nrow(x)
  index = sample(1:n,n_samples,replace = FALSE)
  return(x[index,])
}

list_df_query = lapply(list_df_classes,strat_random_selection, n_samples = 5)

df_query = do.call("rbind",list_df_query)

rm(list_df_query,list_df_classes,legend)

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

read_raster <-function(file_raster, offs, reg){
  geo_info = rgdal::GDALinfo(file_raster,returnStats = FALSE)
  list_bands = list()
  datasource = rgdal::GDAL.open(file_raster)
  for(i in 1:geo_info[3]){
    list_bands[[i]] = rgdal::getRasterData(datasource, offset = offs, region.dim = reg, band=i)
  }
  rgdal::closeDataset(datasource)

  result = array(unlist(list_bands),dim = c(ncols, nrows , geo_info[3]))
  return(result)
}

#example
crop_raster <- function(x,file_raster){
    geo_info = rgdal::GDALinfo(file_raster,returnStats = FALSE)
    coord_pixel = geotransform(geo_info, x)
    xyoffs = coord_pixel$offs
    nrows = coord_pixel$ycount
    ncols = coord_pixel$xcount
    reg = c(nrows,ncols)
    result <- read_raster(file_raster, xyoffs, reg)
    return(result)
}

file_raster = 'IM_20170405_843.tif'

shapes = split(df_query,df_query$OBJECTID)

list_crops = lapply(shapes,crop_raster,file_raster)


a = values(raster('IM_20170405_843.tif'))


#===============================
#Rasterization
#===============================

#rasterize
r <- raster(ncols = reg[1], nrows = reg[2])
extent(r) <- extent(shape)
r <- rasterize(shape, r)


