setwd("~/data/ArcDEM/")
library(sf)
library(tidyverse)
library(dplyr)
library(terra)
library(tidyterra)
library(stars)
library(raster)
library(rgdal) 

############################
# tile
############################
# tile_112 <- rast("29_21_1_1_2m_v3.0_reg_dem.tif") old version
tile_112 <- rast("29_21_1_1_2m_v4.1_dem.tif")
tile_112_mask <- rast("29_21_1_1_2m_v4.1_datamask.tif")
# tile_112_hillshade <- rast("29_21_1_1_2m_v4.1_browse.tif")
plot(tile_112)
plot(tile_112_mask)
# plot(tile_112_hillshade)
# do I need to consider the mask pixel form this mask? 
# tile_122 <- rast("29_21_1_2_2m_v3.0_reg_dem.tif")
tile_122 <- rast("29_21_2_1_2m_v4.1_dem.tif")
tile_122_mask <- rast("29_21_2_1_2m_v4.1_datamask.tif")
plot(tile_122)
plot(tile_122_mask)

tile_DEM <- mosaic(tile_112, tile_122)
plot(tile_DEM)

tile_DEM_proj <- terra::project(tile_DEM, "EPSG:32613")

# tile_crs_string <- crs(tile_DEM, proj=T)

############################
# area of interest
############################
area_interest <- st_read("~/data/areas_of_interest.gpkg")
# area_interest_proj <- st_transform(area_interest, crs=tile_crs_string)
area_interest_proj <- st_transform(area_interest,32613)
plot(tile_DEM_proj)
plot(area_interest_proj, add=T, col="red")

tile_DEM_crop <- crop(tile_DEM_proj, area_interest_proj)
plot(tile_DEM_crop)
plot(area_interest_proj, add=T)

# tile_DEM_crop_proj <- project(tile_DEM_crop, "EPSG:32613")
# area_interest_proj <- st_transform(area_interest, 32613)
# plot(tile_DEM_crop_proj)
# plot(area_interest_proj, add=T)

# writeRaster(tile_DEM_crop, filename = "tile_DEM_crop.tif")
tile_DEM_crop <- rast("~/data/ArcDEM/tile_DEM_crop.tif")
viridis_colors <- viridis::inferno(30)
plot(tile_DEM_crop, col=viridis_colors[30:1])


############################
# Resample
############################
Shannon_map <- rast("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi/SPCA/ALPHA/Shannon_10")
mask <- rast("~/master-thesis/sent_output/mask_sent2_final_NA.tif")
tile_DEM_mask_resample <- terra::project(tile_DEM_crop, mask, method="bilinear")
plot(tile_DEM_mask_resample)
tile_DEM_aggregated <- aggregate(tile_DEM_mask_resample, 10, fun="sd") 
plot(tile_DEM_aggregated)
tile_DEM_aggregated_masked <- mask(tile_DEM_aggregated, Shannon_map)
plot(tile_DEM_aggregated_masked)
plot(Shannon_map, add=T)
table(is.na(tile_DEM_aggregated_masked[1:120540]))
table(is.na(Shannon_map[1:120540]))
Shannon_map_model_DEM <- mask(Shannon_map, tile_DEM_aggregated_masked)
plot(Shannon_map_model_DEM)
plot(tile_DEM_aggregated_masked, add=T, col="red")
table(is.na(tile_DEM_aggregated_masked[1:120540]))
table(is.na(Shannon_map_model_DEM[1:120540]))


Sentinel_map_poly <- as.polygons(Shannon_map_model_DEM, round=F, aggregate=F, extent=F)
Arc_DEM_poly <- as.polygons(tile_DEM_aggregated_masked, round=F, aggregate=F, extent=F)

# writeVector(Sentinel_map_poly, filename = "~/data/output/Sentinel_map_poly.shp")
# writeVector(Arc_DEM_poly, filename = "~/data/output/Arc_DEM_poly.shp")

zero_color <- "red"
# Plot the SpatRaster with data value = 0 in red color
plot(Shannon_map, col = c(zero_color, grey.colors(25)))

Sentinel_vect <- vect("~/data/output/Sentinel_map_poly.shp")
Arc_DEM_vect <- vect("~/data/output/Arc_DEM_poly.shp")
sd_topo <- Arc_DEM_vect$X29_21_1_1
Sentinel_vect$sd_topo <- sd_topo
Sentinel_vect_no_0 <- Sentinel_vect[-c(which(Sentinel_vect$Shannon_10==0))]
spatial <- 1:length(Sentinel_vect_no_0$Shannon_10)
Sentinel_vect_no_0$spatial <- as.numeric(spatial)
# Sentinel_vect_no_0 <- Sentinel_vect[-c(which(Sentinel_vect$Shannon_10==0)),]
# writeVector(Sentinel_vect, filename = "~/data/output/model_object.shp", overwrite=T)
# writeVector(Sentinel_vect_no_0, filename = "~/data/output/model_object_no_0.shp", overwrite=T)
Sentinel_lattice <-readOGR("~/data/output/model_object.shp")


# Sentinel_sf_poly <- st_read("~/data/output/Sentinel_map_poly.shp")
# Arc_DEM_sf_poly <- st_read("~/data/output/Arc_DEM_poly.shp")
# identical(Sentinel_sf_poly$geometry, Arc_DEM_sf_poly$geometry)
# 
# library(rgdal)
# Sentinel_sp_poly <-readOGR("~/data/output/Sentinel_map_poly.shp")
# Arc_DEM_sp_poly <- readOGR("~/data/output/Arc_DEM_poly.shp")
# model_sp <- sp::merge(Sentinel_sp_poly, Arc_DEM_sp_poly)


model_df <- merge(Sentinel_sf_poly, Arc_DEM_sf_poly, by="geometry")
cbind(Sentinel_sf_poly)


#################
# other resampling method
#################
tile_DEM_100m <- aggregate(tile_DEM_crop, fact=100/2, fun="sd") #1.993736
Shannon_map <- rast("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi/SPCA/ALPHA/Shannon_10")
Shannon_proj_DEM <- terra::project(Shannon_map, tile_DEM_100m, method="bilinear")
tile_DEM_100m_masked <- mask(tile_DEM_100m, Shannon_proj_DEM)
plot(Shannon_proj_DEM)
plot(tile_DEM_100m_masked, add=T, col="red")

Shannon_masked <- mask(Shannon_proj_DEM, tile_DEM_100m_masked)
plot(Shannon_masked)
plot(tile_DEM_100m_masked, add=T, col="red")
plot(Shannon_masked$Shannon_10, tile_DEM_100m_masked$`29_21_1_1_2m_v3.0_reg_dem`)

Sentinel_map_poly <- as.polygons(Shannon_masked, round=F, aggregate=F, extent=F)
Arc_DEM_poly <- as.polygons(tile_DEM_100m_masked, round=F, aggregate=F, extent=F)

writeVector(Sentinel_map_poly, filename = "~/data/ArcDEM/Sentinel_map_poly.shp")
writeVector(Arc_DEM_poly, filename = "~/data/ArcDEM/Arc_DEM_poly.shp", overwrite=T)

Sentinel_vect <- vect("~/data/ArcDEM/Sentinel_map_poly.shp")
Arc_DEM_vect <- vect("~/data/ArcDEM/Arc_DEM_poly.shp")
sd_topo <- Arc_DEM_vect$X29_21_1_1
Sentinel_vect$sd_topo <- sd_topo
# Sentinel_vect_no_0 <- Sentinel_vect[-c(which(Sentinel_vect$Shannon_10==0))]
spatial <- 1:length(Sentinel_vect$Shannon_10)
Sentinel_vect$spatial <- as.numeric(spatial)
writeVector(Sentinel_vect, filename = "~/data/ArcDEM/model_object.shp", overwrite=T)
# writeVector(Sentinel_vect_no_0, filename = "~/data/output/model_object_no_0.shp", overwrite=T)
Sentinel_lattice <-readOGR("~/data/ArcDEM/model_object.shp")





