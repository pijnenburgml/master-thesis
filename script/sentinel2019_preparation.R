setwd("~/data")
library(sf)
library(tidyverse)
library(dplyr)
library(terra)
library(tidyterra)
library(stars)
library(raster)
library(cowplot)
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
# read band at 20 m and disaggregate
sent20 <- lapply(1:length(data_list), function(x){rast(data_list[x])})
sent_disagg <- lapply(1:length(sent20), function(x){disagg(sent20[[x]], fact=2, method="bilinear")})

# read band at 10 m resolution
path <- c("/home/mpijne/data/L2A_T13WDS_A012480_20190727T190134_IMG_DATA/R10m", "/home/mpijne/data/L2A_T13WES_A012480_20190727T190134_IMG_DATA/R10m")
data_list <- list.files(path=path, recursive = T, full.names = T, pattern="B.._10m.jp2$")
sent10 <- lapply(1:length(data_list), function(x) {rast(data_list[x])})

# put together
sent10_tot <- c(sent_disagg, sent10)
name <- c()
for(x in 1:length(sent10_tot)){
  name[x] <- names(sent10_tot[[x]])
}
names(sent10_tot) <- name
name <- sort(name)

# mosaic 
sent_aoi <- list()
for (x in 1:(length(name)/2)) {
  pat <- substr(name[x], 24, 26)
  ind <- grep(pattern=pat, names(sent10_tot))
  sent_aoi[[x]] <- mosaic(sent10_tot[[ind[1]]], sent10_tot[[ind[2]]])
}
sent_aoi_stack <- rast(sent_aoi)
plotRGB(sent_aoi_stack, r=3, g=2, b=1, scale=10000, stretch="lin")
writeRaster(sent_aoi_stack, filename="output/full_stack_view.tif")
sent_aoi_stack_crop <- crop(sent_aoi_stack, area_interest_proj)
# writeRaster(sent_aoi_stack_crop, filename="output/sent_crop_view.tif")
# writeRaster(sent_aoi_stack_crop, filetype="ENVI", gdal="INTERLEAVE=BIL", filename = "output/sent_crop_envi", overwrite=T)
sent_aoi_stack_crop <- rast("~/data/output/sent_crop_view.tif")
plotRGB(sent_aoi_stack_crop, r=3, g=2, b=1, scale=10000, stretch="lin")
plot(area_interest_proj, add=T)

################
# data cleaning 
################

# removal of water bodies
names(sent_aoi_stack_crop)
NDWI <- (sent_aoi_stack_crop["T13WDS_20190727T185929_B03_10m"] - sent_aoi_stack_crop["T13WDS_20190727T185929_B08_10m"])/(sent_aoi_stack_crop["T13WDS_20190727T185929_B03_10m"] +sent_aoi_stack_crop["T13WDS_20190727T185929_B08_10m"])
plot(NDWI)
# some water pixel have actually value very close to 0, which is weird but also means that they are not filtered
water_mask <- ifel(NDWI>=0.2, 0, 1)
# writeRaster(water_mask, filename = "~/master-thesis/sent_output/water_mask.tif", overwrite=TRUE)


# shade mask

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


# remove cell below the aoi

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

# Get the intersection coordinate from the two sf objects 
intersc <- st_intersection(sent_poly_sf, area_interest_proj_poly)
text(st_coordinates(intersc)[,1:2], labels=c(1:110))
#from 4 to 97 --> our segment of interest
sep_segment <- matrix(st_coordinates(intersc)[4:97,1:2], ncol=2) # make a matrix form the coordinate of interest
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
mask_below <- mask(sent_aoi_stack_crop, sep_poly_vect, inverse=F, updatevalue=0)
mask_below_try <- ifel(mask_below!=0, 0, 1)
plotRGB(sent_aoi_stack_crop, r=3, g=2, b=1, scale=10000, stretch="lin")
plot(mask_below_try["T13WDS_20190727T185929_B02_10m"], add=T)
# writeRaster(mask_below_try$T13WDS_20190727T185929_B02_10m, filename = "~/data/output/build_are_mask.tif", overwrite=T)


# HAVE TO RUN THE FUNCTION IN THE FILE NAMED HATCHED_FUNCTION
sent_aoi_stack_crop_scaled <- sent_aoi_stack_crop/6
hatched <- patternLayer(x=sep_poly_poly, pattern = "right2left", mode="sfc")
plotRGB(sent_aoi_stack_crop, r=3, g=2, b=1, scale=10000, stretch="lin")
plot(area_interest_proj, add=T, col=sf.colors(categorical=T, alpha=0.6))
plot(hatched, add=T, col="red", lwd=4)

########
# make only 1 mask file out of the several I have
########

mask_wa_sh <- mosaic(water_mask, shade_mask, fun="min")
mask_final <- mosaic(mask_wa_sh, mask_below_try["T13WDS_20190727T185929_B02_10m"], fun="min")
plot(mask_final)
# writeRaster(mask_final, filename = "~/master-thesis/sent_output/mask_sent2_final.tif")

mask <- rast("~/master-thesis/sent_output/mask_sent2_final.tif")
plot(mask)
mask <- ifel(mask==0, NA, 1)
plot(mask)
# writeRaster(mask, filename = "~/data/output/mask_sent2_final_NA.tif")
# writeRaster(mask, filename = "~/data/biodivmapR_sent/mask_sent2_final_NA", filetype="ENVI", gdal="INTERLEAVE=BSQ", overwrite=T, datatype="INT1U")
