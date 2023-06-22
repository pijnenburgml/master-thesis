#!/usr/bin/env R
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=48:00:00

setwd("~/scratch/")
library(sf)
library(tidyverse)
library(tidyterra)
library(dplyr)
library(terra)
library(caTools)
######
# 
######
poly <- vect("~/scratch/polygon_aviris_crop.shp")
poly_proj <- project(poly,"epsg:32613")
poly_ext <- ext(poly_proj)
sent_aoi_stack_crop <- rast("~/data/output/sent_crop_view.tif")
strip_1 <- rast("reflectance_data/ang20190801t160747_rfl")
strip_1_crop <- rectify(strip_1, aoi=poly_ext, filename="stripe_1_crop.tif", overwrite=TRUE)


poly_small <- vect("~/scratch/try2.shp")
poly_proj_small <- project(poly_small,"epsg:32613")
plotRGB(sent_aoi_stack_crop, r=3, g=2, b=1, scale=10000, stretch="lin")
plot(poly_proj_small, col="red", add=T)
strip_crop_small <- crop(strip_1, poly_proj_small)
plot(ext(strip_1), col="red", add=T)
plot(poly_proj_small,  add=T, col="red")

plot(band_1_crop["562.039576 Nanometers"], add=T)
ggplot()+
  geom_spatraster_rgb(data=band_1_crop, interpolate = T)
  

plotRGB(band_1_crop, )

band_1_crop <- crop(band_1, poly_proj)
writeRaster(band_1_crop, filename = "~/scratch/band_1_crop")
plotRGB(sent_aoi_stack_crop, r=3, g=2, b=1, scale=10000, stretch="lin")
plot(ext(band_1), add=T)
plot(poly_proj, add=T, col="red")


