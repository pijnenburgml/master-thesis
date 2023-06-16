setwd("~/data")
library(sf)
library(tidyverse)
library(dplyr)
library(terra)
library(tidyterra)
library(stars)
library(raster)
library(biodivMapR)

############################
# area of interest
############################
area_interest <- st_read("areas_of_interest.gpkg")
area_interest_proj <- st_transform(area_interest,32613)
 
############################
# Sentinel2 data from 11th of August 2019
############################

# need first the spatial downscaling to have all bands at 10m resolution

path <- c("/home/mpijne/data/L2A_T13WDS_A012480_20190727T190134_IMG_DATA/R20m", "/home/mpijne/data/L2A_T13WES_A012480_20190727T190134_IMG_DATA/R20m")
data_list <- list.files(path=path, recursive = T, full.names = T, pattern="B.._20m.jp2$")

# remove bands that are already at 10m resolution

data_list <- data_list[-grep("B02|B03|B04|B08", data_list)]
sent20 <- lapply(1:length(data_list), function(x){rast(data_list[x])})
sent_disagg <- lapply(1:length(sent20), function(x){disagg(sent20[[x]], fact=2, method="bilinear")})
# sent20_crop <- lapply(1:length(sent_disagg), function(x){crop(sent_disagg[[x]], area_interest_proj)})

path <- c("/home/mpijne/data/L2A_T13WDS_A012480_20190727T190134_IMG_DATA/R10m", "/home/mpijne/data/L2A_T13WES_A012480_20190727T190134_IMG_DATA/R10m")
data_list <- list.files(path=path, recursive = T, full.names = T, pattern="B.._10m.jp2$")
sent10 <- lapply(1:length(data_list), function(x) {rast(data_list[x])})
# sent10_crop <- lapply(1:length(sent10), function(x){crop(sent10[[x]], area_interest_proj)})

sent10_tot <- c(sent_disagg, sent10)

name <- c()
for(x in 1:length(sent10_tot)){
  name[x] <- names(sent10_tot[[x]])
}
names(sent10_tot) <- name
name <- sort(name)

sent_aoi <- list()
for (x in 1:(length(name)/2)) {
  pat <- substr(name[x], 24, 26)
  ind <- grep(pattern=pat, names(sent10_tot))
  sent_aoi[[x]] <- mosaic(sent10_tot[[ind[1]]], sent10_tot[[ind[2]]])
}

sent_aoi_stack <- rast(sent_aoi)
plotRGB(sent_aoi_stack, r=3, g=2, b=1, scale=10000, stretch="lin")
lines(area_interest_proj, col="red")
# why lines doesn't actually add a line???
# writeRaster(sent_aoi_stack, filename="output/full_stack_view.tif")

sent_aoi_stack <- rast("output/full_stack_view.tif")
# mkdir ~/master-thesis/sent_output
# writeRaster(sent_aoi_stack, filename = "~/master-thesis/sent_output/sent_aoi_stack.tif", overwrite=T)

sent_aoi_stack_crop <- crop(sent_aoi_stack, area_interest_proj)
# writeRaster(sent_aoi_stack_crop, filename="output/sent_crop_view.tif")
# writeRaster(sent_aoi_stack_crop, filename = "~/master-thesis/sent_output/sent_aoi_stack_crop.tif")

sent_aoi_stack_crop <- rast("output/sent_crop_view.tif")

plotRGB(sent_aoi_stack_crop, r=3, g=2, b=1, scale=10000, stretch="lin")
# why do we have such a different color after cropping 
plot(area_interest_proj, add=T, col="red")


#data cleaning 
#removal of water bodies

names(sent_aoi_stack_crop)
# SWIR <- (sent_aoi_stack_crop["T13WDS_20190727T185929_B11_20m"]+sent_aoi_stack_crop["T13WDS_20190727T185929_B12_20m"])/2
# NDWI<- (SWIR - sent_aoi_stack_crop["T13WDS_20190727T185929_B08_10m"])/(SWIR +sent_aoi_stack_crop["T13WDS_20190727T185929_B08_10m"])
# this version is used for leaf water content

NDWI <- (sent_aoi_stack_crop["T13WDS_20190727T185929_B03_10m"] - sent_aoi_stack_crop["T13WDS_20190727T185929_B08_10m"])/(sent_aoi_stack_crop["T13WDS_20190727T185929_B03_10m"] +sent_aoi_stack_crop["T13WDS_20190727T185929_B08_10m"])
plot(NDWI)
# writeRaster(NDWI, filename = "~/master-thesis/sent_output/NDWI.tif")

# some water pixel have actually value very close to 0, which is weird but also means that they are not filtered
water_mask <- ifel(NDWI>=0.2, 0, 1)
# writeRaster(water_mask, filename = "~/master-thesis/sent_output/water_mask.tif", overwrite=TRUE)

# # see how I can remove the weird pixel in the water that are not removed by the water mask
# mask_test <- ifel(sent_aoi_stack_crop["T13WDS_20190727T185929_B03_10m"]==1 & sent_aoi_stack_crop["T13WDS_20190727T185929_B04_10m"]==1, 1, NA)
# plot(mask_test)
# mask_test_df <- as.data.frame(mask_test, xy=T)
# weird_pixel <- extract(sent_aoi_stack_crop, y=mask_test_df[,1:2])
# mask_test <- ifel(sent_aoi_stack_crop["T13WDS_20190727T185929_B03_10m"]==1 & sent_aoi_stack_crop["T13WDS_20190727T185929_B04_10m"]==1, NA, 1)
# sent_aoi_stack_crop_nowa <- mask(sent_aoi_stack_crop, mask=water_mask)
# sent_aoi_stack_crop_nowa_noweird <- mask(sent_aoi_stack_crop_nowa, mask=mask_test)
# plotRGB(sent_aoi_stack_crop_nowa_noweird, r=3, g=2, b=1
#         ,scale=10000
# )  
# #why so different now?
# # writeRaster(mask_test, filename = "~/master-thesis/sent_output/weird_pixel_mask.tif")
# the weird pixel are removed by the shade mask as well 

# ########
# # NDVI mask?
# ########
# NDVI <- (sent_aoi_stack_crop["T13WDS_20190727T185929_B08_10m"]- sent_aoi_stack_crop["T13WDS_20190727T185929_B04_10m"])/(sent_aoi_stack_crop["T13WDS_20190727T185929_B08_10m"]+sent_aoi_stack_crop["T13WDS_20190727T185929_B04_10m"])
# plot(NDVI)
# NDVI_mask <- ifel(NDVI>-0.5, 1, 0)
# plot(NDVI_mask)
# # writeRaster(NDVI_mask, filename = "~/master-thesis/sent_output/NDVI_mask.tif")
# # The mask doesn't seems to add anything compared to the water + doens't remove the very white looking area like civilized.
# sent_aoi_stack_crop_nowa_NDVI <- mask(sent_aoi_stack_crop_nowa, NDVI_mask)
# ggplot()+
#   geom_spatraster_rgb(data=sent_aoi_stack_crop_nowa_NDVI, r=3, g=2, b=1)
# 
# #######
# # snow pixel removal
# #######
# #removal of snowed pixel
# NDSI <- (sent_aoi_stack_crop_nowa["T13WDS_20190727T185929_B03_10m"] - sent_aoi_stack_crop_nowa["T13WDS_20190727T185929_B11_20m"])/(sent_aoi_stack_crop_nowa["T13WDS_20190727T185929_B03_10m"] + sent_aoi_stack_crop_nowa["T13WDS_20190727T185929_B11_20m"])
# plot(NDSI)
# snow_mask <- ifel(NDSI>=0.4, 0, 1)
# plot(snow_mask)
# writeRaster(snow_mask,filename = "~/master-thesis/sent_output/snow_mask.tif", overwrite=T)


#######
# shade removal
#######

# look at the NIR range of the tundra pixel 
grep(pattern="B08", names(sent_aoi_stack_crop))
for(x in 1:nlyr(sent_aoi_stack_crop)){
  n <- names(sent_aoi_stack_crop)[x]
  if (grepl("B08_10m", n)==T) {
    hist(sent_aoi_stack_crop[[x]]) 
  }
}
shade_mask <- ifel(sent_aoi_stack_crop["T13WDS_20190727T185929_B08_10m"]<1000, 0, 1)
par(mfrow=c(1,2))
plot(shade_mask)
plot(water_mask)
writeRaster(shade_mask, filename = "~/master-thesis/sent_output/shade_mask.tif")


#######
# build area removal
#######

build_area_2022 <- vect("~/master-thesis/sent_output/build_area_polygon.shp")
# ggplot()+
#   geom_spatraster_rgb(data=sent_aoi_stack_crop, r=3, g=2, b=1)+
#   geom_spatvector(data=build_area_2022)

build_pix <- extract(sent_aoi_stack_crop, build_area_2022, xy=T)
# par(mfrow=c(3,3))
# for (x in 2:length(build_pix)) {
#   hist(build_pix[,x], main=paste(x))  
# }
# random_samp <- spatSample(sent_aoi_stack_crop_nowa_noweird, size=16000,method="random", na.rm=T)
# for (x in 1:length(random_samp)) {
#   hist(random_samp[,x], main=paste(x))  
# }

#band 3, 4, 6, 7, 8A, 12 have different distribution, with build area being shifted to higher values.
par(mfrow=c(1,1))
civilized_pix <- ifel(sent_aoi_stack_crop["T13WDS_20190727T185929_B04_10m"]>2000 & sent_aoi_stack_crop["T13WDS_20190727T185929_B03_10m"]>2000 & sent_aoi_stack_crop["T13WDS_20190727T185929_B05_20m"]>2000, 1, NA)
plot(civilized_pix)

sent_aoi_stack_crop[build_area_2022]
try_build_mask <- mask(sent_aoi_stack_crop, build_area_2022, inverse=T, updatevalue=0)
try_build_mask <- ifel(try_build_mask!=0, NA, 0)
plot(try_build_mask)
# ggplot()+
#   geom_spatraster_rgb(data=try_build_mask, r=3, g=2, b=1)
#   # + geom_spatvector(data=build_area_2022)


#######
# remove cell below the aoi
#######

# aoi as an st object
area_interest_proj_poly <- vect(area_interest_proj)
area_interest_proj_poly <- st_as_sf(area_interest_proj_poly)
# the boundaries (rectangle) of the sentinel data as an st object 
sent_poly <- as.polygons(ext(sent_aoi_stack_crop))
sent_poly_sf <- st_as_sf(sent_poly)
st_crs(sent_poly_sf) <- crs(sent_aoi_stack_crop)
plot(sent_poly_sf)
text(st_coordinates(sent_poly_sf)[,1:2], labels=c(1:5))
# 476480, 7657730; 517470, 7657730 --> coordinates of the two bottom edge 

# Get the intersection coordinate from the tow sf object 
try3 <- st_intersection(sent_poly_sf, area_interest_proj_poly)
text(st_coordinates(try3)[,1:2], labels=c(1:110))
#from 4 to 97 --> our segment of interest
sep_segment <- matrix(st_coordinates(try3)[4:97,1:2], ncol=2) # make a matrix form the coordinate of interest
coin <- matrix(c(476480, 7657730, 515773.2, 7657730), ncol=2, byrow = T) # make a matrix with the coordinate of the bottom left edge + the first coordinate of the segment to be able to form a closed polygon 
sep_segment <- rbind(sep_segment, coin) #bind together 
plot(sep_segment)
sep_poly <- st_linestring(sep_segment) # make a linestring object out of our matrix
sep_poly_poly <- st_cast(sep_poly, "POLYGON") # make a polygon out of the linestring
plot(sep_poly_poly) 
sep_poly_vect <- vect(sep_poly_poly) # make a SpatVector object
crs(sep_poly_vect) <- crs(sent_aoi_stack_crop) 
plotRGB(sent_aoi_stack_crop, r=3, g=2, b=1, scale=10000, stretch="lin")
plot(sep_poly_vect, add=T, col="green")

mask_below <- mask(sent_aoi_stack_crop, sep_poly_vect, inverse=T, updatevalue=0)
plot(mask_below)
mask_below_try <- ifel(mask_below!=0, NA, 0)
plotRGB(sent_aoi_stack_crop, r=3, g=2, b=1, scale=10000, stretch="lin")
plot(mask_below_try["T13WDS_20190727T185929_B02_10m"], add=T)

########
# make only 1 mask file out of the several I have
########

mask_wa_sh <- mosaic(water_mask, shade_mask, fun="min")
mask_wa_sh_bul <- mosaic(mask_wa_sh, try_build_mask["T13WDS_20190727T185929_B02_10m"], fun="min")
plot(mask_wa_sh_bul)

mask_final <- mosaic(mask_wa_sh_bul, mask_below_try["T13WDS_20190727T185929_B02_10m"], fun="min")
plot(mask_final)
# writeRaster(mask_final, filename = "~/master-thesis/sent_output/mask_sent2_final.tif")


