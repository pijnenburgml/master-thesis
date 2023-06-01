setwd("C:/Users/Marie/OneDrive/Documents/master/master thesis/master-thesis")
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
# Sentinel2 data from 11th of August 2017
############################

# #Mosaic 2 tiles together:
path <- c("C:/Users/Marie/OneDrive/Documents/master/master thesis/master-thesis/data 2017/S2A_MSIL1C_20170811T185921_N0205_R013_T13WES_20170811T190305.SAFE", "C:/Users/Marie/OneDrive/Documents/master/master thesis/master-thesis/data 2017/S2A_MSIL1C_20170811T185921_N0205_R013_T13WDS_20170811T190305.SAFE")
#WDT outside of the area of interest
data_list <- list.files(path=path, recursive = T, full.names = T, pattern = "B...jp2$")
data_list <- data_list[-grep("B01|B09|B10", data_list)]
sent <- lapply(1:length(data_list), function(x) {rast(data_list[x])})
sent_tot <- list()
for (x in 1:(length(sent)/2)) {
    sent_tot[x] <- mosaic(sent[[x]], sent[[x+length(sent)/2]])
}
sent_tot
#Aggregating 10m resulotion band to 20 meter resolution bands:
ten_band <- c(1,2,3,7)
sent_aggr <- lapply(ten_band, function(x) {aggregate(sent_tot[[x]], 2)})
sent_20 <- c(sent_aggr[c(1,2,3)], sent_tot[c(4,5,6)], sent_aggr[4], sent_tot[c(8,9,10)])
sent_20_stack <- rast(sent_20)
sent_20_stack
ext(area_interest_proj)
sent_20_stack_int <- crop(sent_20_stack, area_interest_proj)
writeRaster(sent_20_stack_int, "sentinel2_aoi_20mr.raw")
# #writeRaster(sent_20_stack_int, "sentinel2_aoi_20mr.tif")
sentinel2_aoi_20mr <- rast("sentinel2_aoi_20mr.tif")
summary(sentinel2_aoi_20mr)
plotRGB(sentinel2_aoi_20mr, r=3, g=2, b=1 ,scale=10000, stretch="lin")
lines(area_interest_proj, col="red")


# plotRGB(sent_20_stack_int, r=3, g=2, b=1 ,scale=10000, stretch="lin")
# lines(area_interest_proj, col="red")
# #crop worked well


#data cleaning 
#removal of water bodies

names(sent_20_stack_int)
SWIR <- (sent_20_stack_int["T13WDS_20170811T185921_B11"]+sent_20_stack_int["T13WDS_20170811T185921_B12"])/2
#NDBI <- (SWIR - sent_20_stack_int["T13WDS_20170811T185921_B08"])/(SWIR +sent_20_stack_int["T13WDS_20170811T185921_B08"])
NDBI <- (sent_20_stack_int["T13WDS_20170811T185921_B03"] - sent_20_stack_int["T13WDS_20170811T185921_B08"])/(sent_20_stack_int["T13WDS_20170811T185921_B03"] +sent_20_stack_int["T13WDS_20170811T185921_B08"])
water_mask <- ifel(NDBI>=0.2, NA, 1)
plot(NDBI)
sent_20_stack_int_nowa <- mask(sent_20_stack_int, mask=water_mask)
plotRGB(sent_20_stack_int_nowa, r=3, g=2, b=1 ,scale=10000, stretch="lin")  
#why so different now?

#removal of snowed pixel
NDSI <- (sent_20_stack_int_nowa["T13WDS_20170811T185921_B03"] - sent_20_stack_int_nowa["T13WDS_20170811T185921_B11"])/(sent_20_stack_int_nowa["T13WDS_20170811T185921_B03"] + sent_20_stack_int_nowa["T13WDS_20170811T185921_B11"])
plot(NDSI)

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




