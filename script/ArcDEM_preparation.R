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
# writeRaster(tile_DEM_proj, "~/data/ArcDEM/tile_DEM_proj.tif")
# tile_crs_string <- crs(tile_DEM, proj=T)

plot(tile_DEM_proj, colNA="red")


############################
# area of interest
############################
# area_interest <- st_read("~/data/areas_of_interest.gpkg")
# # area_interest_proj <- st_transform(area_interest, crs=tile_crs_string)
# area_interest_proj <- st_transform(area_interest,32613)
# plot(tile_DEM_proj)
# plot(area_interest_proj, add=T, col="red")
# 
# tile_DEM_crop <- crop(tile_DEM_proj, area_interest_proj)
# plot(tile_DEM_crop)
# plot(area_interest_proj, add=T)
# 
# # tile_DEM_crop_proj <- project(tile_DEM_crop, "EPSG:32613")
# # area_interest_proj <- st_transform(area_interest, 32613)
# # plot(tile_DEM_crop_proj)
# # plot(area_interest_proj, add=T)
# 
# # writeRaster(tile_DEM_crop, filename = "tile_DEM_crop_v4.tif")
# tile_DEM_crop <- rast("~/data/ArcDEM/tile_DEM_crop_v4.tif")
# viridis_colors <- viridis::inferno(30)
# plot(tile_DEM_crop, col=viridis_colors[30:1])
# plot(tile_DEM_crop, colNA="red")

#################
# try to crop larger before aggregating
#################

area_interest <- st_read("~/data/areas_of_interest.gpkg")
area_interest_proj <- st_transform(area_interest,32613)
ext(area_interest_proj)
extended_aoi <- ext(ext(area_interest_proj)[1]-100, ext(area_interest_proj)[2]+100, ext(area_interest_proj)[3]-100, ext(area_interest_proj)[4]+100)
crs(area_interest_proj)
tile_DEM_crop <- crop(tile_DEM_proj, extended_aoi)
plot(tile_DEM_crop)
ext(tile_DEM_proj)
ext(area_interest_proj)
Shannon_map <- rast("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi/SPCA/ALPHA/Shannon_10")
temp_rast <- rast(ext(Shannon_map), resolution = 2)
tile_DEM_crop_resample <- resample(tile_DEM_crop, temp_rast, method = "bilinear")
plot(tile_DEM_crop_resample, colNA="red")
tile_DEM_crop_resample_noNA <- subst(tile_DEM_crop_resample, NA, -33)
tile_DEM_crop_aggregated <- aggregate(tile_DEM_crop_resample_noNA, fact=100/2, fun="sd")
plot(tile_DEM_crop_aggregated)
tile_DEM_crop_masked <- mask(tile_DEM_crop_aggregated, Shannon_map)
plot(tile_DEM_crop_masked, col="red")
plot(Shannon_map, add=T)
ncell(Shannon_map)
ncell(tile_DEM_crop_masked)

# produce object for the model
table(is.na(tile_DEM_crop_masked[1:120540]))
table(is.na(Shannon_map[1:120540]))

Sentinel_map_poly <- as.polygons(Shannon_map, round=F, aggregate=F, extent=F)
Arc_DEM_poly <- as.polygons(tile_DEM_crop_masked, round=F, aggregate=F, extent=F)
# writeVector(Sentinel_map_poly, filename = "~/data/output/Sentinel_map_poly.shp", overwrite=T)
# writeVector(Arc_DEM_poly, filename = "~/data/output/Arc_DEM_poly.shp", overwrite=T)
Sentinel_vect <- vect("~/data/output/Sentinel_map_poly.shp")
Arc_DEM_vect <- vect("~/data/output/Arc_DEM_poly.shp")
sd_topo <- Arc_DEM_vect$X29_21_1_1
Sentinel_vect$sd_topo <- sd_topo
# Sentinel_vect_no_0 <- Sentinel_vect[-c(which(Sentinel_vect$Shannon_10==0))]
# Sentinel_vect_no_0 <- Sentinel_vect[-c(which(Sentinel_vect$Shannon_10==0)),]

# writeVector(Sentinel_vect, filename = "~/data/output/model_object.shp", overwrite=T)
# writeVector(Sentinel_vect_no_0, filename = "~/data/output/model_object_no_0.shp", overwrite=T)
Sentinel_lattice <-readOGR("~/data/output/model_object.shp")

############################
# Resample (best way)
############################
tile_DEM_crop <- rast("~/data/ArcDEM/tile_DEM_crop_v4.tif")
Shannon_map <- rast("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi/SPCA/ALPHA/Shannon_10")
plot(Shannon_map, alpha = 0.5, type="classes")
plot(tile_DEM_crop, alpha = 0.5, add = T, type="classes")
temp_rast <- rast(ext(Shannon_map), resolution = 2)
tile_DEM_crop_resample <- resample(tile_DEM_crop, temp_rast, method = "bilinear")
plot(tile_DEM_crop_resample, colNA="red")
# no NA at the border here

ext(tile_DEM_crop_resample)[c(1,3)]+5
plot(tile_DEM_crop_resample)
terra::points(ext(tile_DEM_crop_resample)[1], col="red", cex=2, add=T)
tile_DEM_crop_resample_noNA <- subst(tile_DEM_crop_resample, NA, -32)
plot(tile_DEM_crop_resample_noNA, colNA="red")

# point_coords <- origin(tile_DEM_crop_resample)
# point <- vect(st_point(point_coords))
# plot(tile_DEM_crop_resample, colNA="red")
# points(point, col="blue", cex=10)

# plot(Shannon_map, add=T, col="grey", colNA="red")
# point_coords <- origin(Shannon_map)
# point <- vect(st_point(point_coords))
# plot(Shannon_map)
# points(point, col="black", cex=7)

tile_DEM_crop_aggregated <- aggregate(tile_DEM_crop_resample, fact=100/2, fun="sd")
tile_DEM_crop_aggregated_try <- aggregate(tile_DEM_crop_resample_noNA, fact=100/2, fun="sd")
plot(tile_DEM_crop_aggregated_try, colNA="red")

ncell(tile_DEM_crop_aggregated)
table(is.na(tile_DEM_crop_aggregated[1:120540]))
plot(tile_DEM_crop_aggregated, colNA="red")
# aggregate function produce NA at the border. 

# option to mask the Shannon to remove the line of NAs 
tile_DEM_crop_masked <- mask(tile_DEM_crop_aggregated, Shannon_map)
plot(Shannon_map)
plot(tile_DEM_crop_masked, add=T, col="red")
Shannon_map_masked <- mask(Shannon_map, tile_DEM_crop_masked)
table(is.na(Shannon_map_masked[1:120540]))
table(is.na(tile_DEM_crop_masked[1:120540]))


#######################
# first way to resample (worse, loose a lot of data)
#######################

tile_DEM_crop <- rast("~/data/ArcDEM/tile_DEM_crop_v4.tif")
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

# produce object for the model
Sentinel_map_poly <- as.polygons(Shannon_map_model_DEM, round=F, aggregate=F, extent=F)
Arc_DEM_poly <- as.polygons(tile_DEM_aggregated_masked, round=F, aggregate=F, extent=F)
# writeVector(Sentinel_map_poly, filename = "~/data/output/Sentinel_map_poly.shp")
# writeVector(Arc_DEM_poly, filename = "~/data/output/Arc_DEM_poly.shp")
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

model_df <- merge(Sentinel_sf_poly, Arc_DEM_sf_poly, by="geometry")
cbind(Sentinel_sf_poly)

# observe where 0 values are
zero_color <- "red"
# Plot the SpatRaster with data value = 0 in red color
plot(Shannon_map, col = c(zero_color, grey.colors(25)))

#########################
# other resampling method (a bit better)
#########################
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

# make an object for modelling
Sentinel_map_poly <- as.polygons(Shannon_masked, round=F, aggregate=F, extent=F)
Arc_DEM_poly <- as.polygons(tile_DEM_100m_masked, round=F, aggregate=F, extent=F)
# writeVector(Sentinel_map_poly, filename = "~/data/ArcDEM/Sentinel_map_poly.shp")
# writeVector(Arc_DEM_poly, filename = "~/data/ArcDEM/Arc_DEM_poly.shp", overwrite=T)

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


#####################
# variogram
#####################
library(gstat)
library(sp)
library(raster)
sd_ele_100 <- rast("~/data/ArcDEM/ArcDEM_masked_10_res.tif")
sd_ele_100_raster <- raster(sd_ele_100)
sd_ele_100_sp <- as(sd_ele_100_raster, "SpatialPointsDataFrame")

vario_ele <- variogram(X29_21_1_1_2m_v4.1_dem~1, data=sd_ele_100_sp, cutoff=30000, width=500)
plot(vario_ele)

vario_ele_close <- variogram(X29_21_1_1_2m_v4.1_dem~1, data=sd_ele_100_sp, cutoff=5000, width=100)
plot(vario_ele_close)

vario_ele_far <- variogram(X29_21_1_1_2m_v4.1_dem~1, data=sd_ele_100_sp, cutoff=30000, width=1000)
plot(vario_ele_far)

vario_log_ele <- variogram(log(X29_21_1_1_2m_v4.1_dem)~1, data=sd_ele_100_sp)
plot(vario_log_ele)

vario_log_ele_close <- variogram(log(X29_21_1_1_2m_v4.1_dem)~1, data=sd_ele_100_sp, cutoff=5000, width=100)
plot(vario_log_ele_close)

vario_log_ele_far <- variogram(log(X29_21_1_1_2m_v4.1_dem)~1, data=sd_ele_100_sp, cutoff=30000, width=1000)
plot(vario_log_ele_far)

try <- vgm(0.1, "Exp", 5000, nugget=0.1)
plot(try, cutoff=15000)

vario_fit_ele <- fit.variogram(vario_log_ele, vgm(0.1, "Exp", 5000, nugget=0.1), fit.kappa = TRUE)
plot(vario_fit_ele, cutoff=10000)
plot(vario_log_ele, vario_fit_ele)
plot(vario_log_ele_far, vario_fit_ele)


vario_fit_close <- fit.variogram(vario_close, vgm(0.05, "Mat", 4000, nugget=0.12), fit.kappa = TRUE)
vario_fit_far <- fit.variogram(vario_far, vgm(0.05, "Mat", 4000, nugget=0.12), fit.kappa = TRUE)

plot(vario_fit_far, cutoff=5000)
plot(vario_far, vario_fit_far)


plot(vario_fit_close, cutoff=5000)
plot(vario_close, vario_fit_close)


