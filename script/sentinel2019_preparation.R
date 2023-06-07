setwd("data")
library(sf)
library(tidyverse)
library(dplyr)
library(terra)
# library(biodivMapR)

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

sent_aoi_stack_crop <- crop(sent_aoi_stack, area_interest_proj)
writeRaster(sent_aoi_stack_crop, filename="output/sent_crop_view.tif")

plotRGB(sent_aoi_stack_crop, r=3, g=2, b=1, scale=10000, stretch="lin")
# why do we have such a different color after cropping 
plot(area_interest_proj, add=T, col="red")


#data cleaning 
#removal of water bodies

names(sent_aoi_stack_crop)
SWIR <- (sent_20_stack_int["T13WDS_20170811T185921_B11"]+sent_20_stack_int["T13WDS_20170811T185921_B12"])/2
#NDWI <- (SWIR - sent_20_stack_int["T13WDS_20170811T185921_B08"])/(SWIR +sent_20_stack_int["T13WDS_20170811T185921_B08"])
NDWI <- (sent_aoi_stack_crop["T13WDS_20190727T185929_B03_10m"] - sent_aoi_stack_crop["T13WDS_20190727T185929_B08_10m"])/(sent_aoi_stack_crop["T13WDS_20190727T185929_B03_10m"] +sent_aoi_stack_crop["T13WDS_20190727T185929_B08_10m"])
plot(NDWI)
# some water pixel have actually value very close to 0, which is weird but also means that they are not filtered
water_mask <- ifel(NDWI>=0.2, NA, 1)
# writeRaster(water_mask, filename = "output/water_mask.tif", overwrite=TRUE)
mask_test <- ifel(sent_aoi_stack_crop["T13WDS_20190727T185929_B03_10m"]==1 & sent_aoi_stack_crop["T13WDS_20190727T185929_B04_10m"]==1, NA, 1)
plot(mask_test)
sent_aoi_stack_crop_nowa <- mask(sent_aoi_stack_crop, mask=water_mask)
plotRGB(sent_aoi_stack_crop_nowa, r=3, g=2, b=1
        ,scale=10000, stretch="lin"
)  
#why so different now?

#removal of snowed pixel
NDSI <- (sent_aoi_stack_crop_nowa["T13WDS_20190727T185929_B03_10m"] - sent_aoi_stack_crop_nowa["T13WDS_20190727T185929_B11_20m"])/(sent_aoi_stack_crop_nowa["T13WDS_20190727T185929_B03_10m"] + sent_aoi_stack_crop_nowa["T13WDS_20190727T185929_B11_20m"])
plot(NDSI)

# removal of shade
# look at the NIR range of the tundra pixel 
grep(pattern="B08", names(sent_aoi_stack_crop))
for(x in 1:nlyr(sent_aoi_stack_crop)){
  n <- names(sent_aoi_stack_crop)[x]
  if (grepl("B08_10m", n)==T) {
    hist(sent_aoi_stack_crop[[x]]) 
  }
}

shade_mask <- ifel(sent_aoi_stack_crop["T13WDS_20190727T185929_B08_10m"]<1000, NA, 1)
par(mfrow=c(1,2))
plot(shade_mask)
plot(water_mask)

########
NDVI <- (sent_aoi_stack_crop["T13WDS_20190727T185929_B08_10m"]- sent_aoi_stack_crop["T13WDS_20190727T185929_B04_10m"])/(sent_aoi_stack_crop["T13WDS_20190727T185929_B08_10m"]+sent_aoi_stack_crop["T13WDS_20190727T185929_B04_10m"])
plot(NDVI)
NDVI_mask <- ifel(NDVI>-0.5, 1, NA)


################################
# Sentinel2 data from 27th of July 2019
################################
############################

#Mosaic 2 tiles together:
path <- c("C:/Users/Marie/OneDrive/Documents/master/master thesis/master-thesis/data 2019/S2B_MSIL1C_20190727T185929_N0208_R013_T13WDS_20190727T210028.SAFE", "C:/Users/Marie/OneDrive/Documents/master/master thesis/master-thesis/data 2019/S2B_MSIL2A_20190727T185929_N0213_R013_T13WES_20190727T214238.SAFE/GRANULE/L2A_T13WES_A012480_20190727T190134/IMG_DATA/R10m")
data_list <- list.files(path=path, recursive = T, full.names = T, pattern = "B*jp2$")
# data_list <- data_list[-grep("B01|B09|B10|TCI|PVI|AOT|WVP", data_list)]
data_list <- data_list[grep("B02|B03|B04|B08", data_list)]
sent <- lapply(1:length(data_list), function(x) {rast(data_list[x])})
sent_tot <- list()
for (x in 1:(length(sent)/2)) {
  sent_tot[x] <- mosaic(sent[[x]], sent[[x+length(sent)/2]])
}
sent_tot
#Aggregating 10m resulotion band to 20 meter resolution bands:
# ten_band <- c(1,2,3,7)
# sent_aggr <- lapply(ten_band, function(x) {aggregate(sent_tot[[x]], 2)})
# sent_20 <- c(sent_aggr[c(1,2,3)], sent_tot[c(4,5,6)], sent_aggr[4], sent_tot[c(8,9,10)])
# sent_20_stack <- rast(sent_20)
# sent_20_stack
# ext(area_interest_proj)
# sent_20_stack_int <- crop(sent_20_stack, area_interest_proj)
#writeRaster(sent_20_stack_int, "sentinel2_aoi_20mr.tif")
sentinel2_aoi_20mr <- rast("sentinel2_aoi_20mr.tif")
summary(sentinel2_aoi_20mr)
plotRGB(sentinel2_aoi_20mr, r=3, g=2, b=1 ,scale=10000, stretch="lin")
lines(area_interest_proj, col="red")
sent_stack <- rast(sent_tot)
plotRGB(sent_stack, r=3, g=2, b=1 ,scale=10000, stretch="lin")
lines(area_interest_proj, col="red")