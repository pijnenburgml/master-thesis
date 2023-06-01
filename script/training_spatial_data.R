setwd("C:/Users/Marie/OneDrive/Documents/master/master thesis/master-thesis")
library(sf)
library(tidyverse)
library(dplyr)
library(terra)
###
#download of the area of interst

area_interest <- st_read("areas_of_interest.gpkg")
str(area_interest)
summary(area_interest)
class(area_interest)
names(area_interest)
plot(area_interest)
#area_interest_tidy <- read_sf("areas_of_interest.gpkg") 
#class(area_interest_tidy)
#summary(area_interest_tidy)
st_crs(area_interest)
crs(area_interest)
area_interest_proj <- st_transform(area_interest,32613)
crs(area_interest_proj)
plot(area_interest_proj)

###
#download sentinel2 data from computer
where <- c("C:/Users/Marie/OneDrive/Documents/master/master thesis/master-thesis/S2A_MSIL2A_20220904T185931_N0400_R013_T13WES_20220905T003602.SAFE")
data_list <- list.files(path=where, recursive = T, full.names = T, pattern = "B.._20m.jp2$")
layers_pic1 <- lapply(1:length(data_list), function(x) {rast(data_list[x])})
stack_layers_pic1 <- rast(layers_pic1)
class(stack_layers_pic1)
plot(stack_layers_pic1)
nlyr(stack_layers_pic1)
#stack_layers_pic1_scaled <- stack_layers_pic1/10000
plotRGB(stack_layers_pic1, r=4, g=3, b=2, scale=10000, stretch="lin")

#second cell
where <- c("C:/Users/Marie/OneDrive/Documents/master/master thesis/master-thesis/S2A_MSIL2A_20220904T185931_N0400_R013_T13WDS_20220905T003602.SAFE")
data_list <- list.files(path=where, recursive = T, full.names = T, pattern = "B.._20m.jp2$")
layers_pic2 <- lapply(1:length(data_list), function(x) {rast(data_list[x])})
stack_layers_pic2 <- rast(layers_pic2)
plotRGB(stack_layers_pic2, r=4, g=3, b=2, scale=10000, stretch="lin")

#use true color image directly
m1 <- rast("C:/Users/Marie/OneDrive/Documents/master/master thesis/master-thesis/S2A_MSIL2A_20220904T185931_N0400_R013_T13WES_20220905T003602.SAFE/GRANULE/L2A_T13WES_A037619_20220904T185928/IMG_DATA/R20m/T13WES_20220904T185931_TCI_20m.jp2")
plot(m1)
m1
#coord. ref. : WGS 84 / UTM zone 13N (EPSG:32613) => to use to reproject area of interest data
crs(m1)
plot(m1)
plot(area_interest_proj, add=T, col="red")

m2 <- rast("C:/Users/Marie/OneDrive/Documents/master/master thesis/master-thesis/S2A_MSIL2A_20220904T185931_N0400_R013_T13WDS_20220905T003602.SAFE/GRANULE/L2A_T13WDS_A037619_20220904T185928/IMG_DATA/R20m/T13WDS_20220904T185931_TCI_20m.jp2")
plot(m2)
m2

mtot <- merge(m1, m2)
plot(mtot)
lines(area_interest_proj, col="red")

mmosaic <- mosaic(m1, m2)
plotRGB(mmosaic)
####
#map a backgroud tile with the area of interest
library(maptiles)
bg <- get_tiles(mtot)
plotRGB(bg)
plot(mtot, add=T)
lines(area_interest_proj, col="red")

#to be continued
# try <- crop(mtot, area_interest_proj)
# plot(try, add=T)
