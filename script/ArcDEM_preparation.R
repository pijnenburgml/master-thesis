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
tile_112 <- rast("29_21_1_1_2m_v4.1_dem.tif")
tile_112_mask <- rast("29_21_1_1_2m_v4.1_datamask.tif")
plot(tile_112)
plot(tile_112_mask)
tile_122 <- rast("29_21_2_1_2m_v4.1_dem.tif")
tile_122_mask <- rast("29_21_2_1_2m_v4.1_datamask.tif")
plot(tile_122)
plot(tile_122_mask)
tile_DEM <- mosaic(tile_112, tile_122)
plot(tile_DEM)
tile_DEM_proj <- terra::project(tile_DEM, "EPSG:32613")
# writeRaster(tile_DEM_proj, "~/data/ArcDEM/tile_DEM_proj.tif")
tile_DEM_proj <- rast("~/data/ArcDEM/tile_DEM_proj.tif")

# production of topography data done directly in the scripts for modelling  

# #################
# # Production of the sd(elevation) data
# #################
# 
# area_interest <- st_read("~/data/areas_of_interest.gpkg")
# area_interest_proj <- st_transform(area_interest,32613)
# ext(area_interest_proj)
# extended_aoi <- ext(ext(area_interest_proj)[1]-100, ext(area_interest_proj)[2]+100, ext(area_interest_proj)[3]-100, ext(area_interest_proj)[4]+100)
# crs(area_interest_proj)
# tile_DEM_crop <- crop(tile_DEM_proj, extended_aoi)
# plot(tile_DEM_crop)
# ext(tile_DEM_proj)
# ext(area_interest_proj)
# Shannon_map <- rast("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi/SPCA/ALPHA/Shannon_10")
# temp_rast <- rast(ext(Shannon_map), resolution = 2)
# tile_DEM_crop_resample <- resample(tile_DEM_crop, temp_rast, method = "bilinear")
# plot(tile_DEM_crop_resample, colNA="red")
# tile_DEM_crop_resample_noNA <- subst(tile_DEM_crop_resample, NA, -33)
# tile_DEM_crop_aggregated <- aggregate(tile_DEM_crop_resample_noNA, fact=100/2, fun="sd")
# plot(tile_DEM_crop_aggregated)
# tile_DEM_crop_masked <- mask(tile_DEM_crop_aggregated, Shannon_map)
# plot(tile_DEM_crop_masked, col="red")
# plot(Shannon_map, add=T)
# ncell(Shannon_map)
# ncell(tile_DEM_crop_masked)
# 
# 
# #####################
# # semi-variogram analysis
# #####################
# library(gstat)
# library(sp)
# library(raster)
# # sd(elevation)
# sd_ele_100 <- rast("~/data/ArcDEM/ArcDEM_masked_10_res.tif")
# sd_ele_100_raster <- raster(sd_ele_100)
# sd_ele_100_sp <- as(sd_ele_100_raster, "SpatialPointsDataFrame")
# vario_log_ele_close <- variogram(log(X29_21_1_1_2m_v4.1_dem)~1, data=sd_ele_100_sp, cutoff=5000, width=100)
# plot(vario_log_ele_close)
# try <- vgm(0.1, "Exp", 1000, nugget=0.12)
# plot(try, cutoff=5000)
# vario_fit_close <- fit.variogram(vario_log_ele_close, try, fit.kappa = TRUE)
# plot(vario_log_ele_close, vario_fit_close, xlab="distance [m]")
# vario_fit_close
# plot(vario_log_ele_close, vario_fit_close)
# preds = variogramLine(vario_fit_close, maxdist = 5000)
# semivar_ele <- ggplot()+
#   geom_point(data=vario_log_ele_close, aes(x=dist, y=gamma), pch=21, cex=2,  col="dodgerblue1")+
#   geom_line(data=preds, aes(x=dist, y=gamma), col="dodgerblue1")+
#   labs(y=expression(paste("semi-variance of ", sigma, ("elevation"))), x=("distance [m]"))+
#   theme_classic()+
#   theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), text=element_text(size=14))+
#   coord_cartesian(xlim =c(250, 5000), ylim = c(0, 0.7))
# 
# # sd slope
# sd_slope_100 <- rast("~/data/ArcDEM/sd_slope_masked_10_res.tif")
# sd_slope_100_raster <- raster(sd_slope_100)
# sd_slope_100_sp <- as(sd_slope_100_raster, "SpatialPointsDataFrame")
# vario_slope_close <- variogram(log(slope)~1, data=sd_slope_100_sp, cutoff=5000, width=100)
# plot(vario_slope_close)
# try <- vgm(0.1, "Exp", 1000, nugget=0.12)
# plot(try, cutoff=5000)
# vario_fit_slope_close <- fit.variogram(vario_slope_close, try, fit.kappa = TRUE)
# plot(vario_slope_close, vario_fit_slope_close, xlab="distance [m]")
# vario_fit_slope_close
# preds = variogramLine(vario_fit_slope_close, maxdist = 5000)
# semivar_slope <- ggplot()+
#   geom_point(data=vario_slope_close, aes(x=dist, y=gamma), pch=21, cex=2,  col="dodgerblue1")+
#   geom_line(data=preds, aes(x=dist, y=gamma), col="dodgerblue1")+
#   labs(y=expression(paste("semi-variance of ", sigma, ("slope"))), x=("distance [m]"))+
#   theme_classic()+
#   theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), text=element_text(size=14))+
#   coord_cartesian(xlim =c(250, 5000), ylim = c(0, 0.4))
# 
# 
# ##################
# # Terrain
# ##################
# tile_DEM_proj <- rast("~/data/ArcDEM/tile_DEM_proj.tif")
# area_interest <- st_read("~/data/areas_of_interest.gpkg")
# area_interest_proj <- st_transform(area_interest,32613)
# ext(area_interest_proj)
# window_size <- c(10, 20, 30, 50, 100)
# multiple <- window_size*10
# extended_aoi <- ext(ext(area_interest_proj)[1]-multiple[5], ext(area_interest_proj)[2]+multiple[5], ext(area_interest_proj)[3]-multiple[5], ext(area_interest_proj)[4]+multiple[5])
# tile_DEM_crop <- crop(tile_DEM_proj, extended_aoi)
# fact <- window_size*10/2
# 
# for(x in 1:length(fact)){
#   # browser()
#   Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", window_size[x], "PC1278", sep="_"))
#   temp_rast <- rast(ext(Shannon), resolution = 2)
#   tile_DEM_crop_resample <- resample(tile_DEM_crop, temp_rast, method = "bilinear")
#   slope <- terrain(tile_DEM_crop_resample, v="slope", neighbors=4, unit="degrees")
#   # aspect <- terrain(tile_DEM_crop_resample, v="aspect", neighbors=4, unit="degrees")
#   slope_noNA <- subst(slope, NA, 0) #have to adjust
#   # aspect_noNA <- subst(aspect, NA, 0)
#   # slope_aggregated <- aggregate(slope_noNA, fact=window_size[x], fun="sd")
#   slope_aggregated <- aggregate(slope_noNA, fact=window_size[x], fun="mean")
#   slope_mask <- mask(slope_aggregated, Shannon)
#   # aspect_mask <- mask(aspect_aggregated, Shannon)
#   writeRaster(slope_mask, filename = paste("~/data/ArcDEM/mean_slope_masked_", window_size[x], "_res.tif", sep=""))
#   # writeRaster(aspect_mask, filename = paste("~/data/ArcDEM/mean_aspect_masked_", window_size[x], "_res.tif", sep=""))
#   # Arc_DEM_poly <- as.polygons(tile_DEM_masked, round=F, aggregate=F, extent=F, na.rm=F)
#   # Sentinel_shannon_poly <- as.polygons(Shannon, round=F, aggregate=F, extent=F,na.rm=F)
#   # writeVector(Sentinel_shannon_poly,filename=paste("~/data/output/INLA_modelling/Sentinel_shannon_poly", "res", window_size[x], "with_NA", sep="_"))
#   # writeVector(Arc_DEM_poly, filename = paste("~/data/output/INLA_modelling/Arc_DEM_poly","res",window_size[x], "with_NA", sep="_"))
# }
# 
# ###############
# # Effect of max elevation
# ###############
# 
# tile_DEM_proj <- rast("~/data/ArcDEM/tile_DEM_proj.tif")
# area_interest <- st_read("~/data/areas_of_interest.gpkg")
# area_interest_proj <- st_transform(area_interest,32613)
# ext(area_interest_proj)
# window_size <- c(10, 20, 30, 50, 100)
# multiple <- window_size*10
# extended_aoi <- ext(ext(area_interest_proj)[1]-multiple[5], ext(area_interest_proj)[2]+multiple[5], ext(area_interest_proj)[3]-multiple[5], ext(area_interest_proj)[4]+multiple[5])
# tile_DEM_crop <- crop(tile_DEM_proj, extended_aoi)
# fact <- window_size*10/2
# 
# for(x in 1:length(fact)){
#   # browser()
#   Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", window_size[x], "PC1278", sep="_"))
#   temp_rast <- rast(ext(Shannon), resolution = 2)
#   tile_DEM_crop_resample <- resample(tile_DEM_crop, temp_rast, method = "bilinear")
#   tile_DEM_crop_resample_noNA <- subst(tile_DEM_crop_resample, NA, -33) #have to adjust
#   tile_DEM_crop_aggregated <- aggregate(tile_DEM_crop_resample_noNA, fact=fact[x], fun="max")
#   tile_DEM_masked <- mask(tile_DEM_crop_aggregated, Shannon)
#   writeRaster(tile_DEM_masked, filename = paste("~/data/ArcDEM/max_ele_masked_", window_size[x], "_res.tif", sep=""))
#   # Arc_DEM_poly <- as.polygons(tile_DEM_masked, round=F, aggregate=F, extent=F, na.rm=F)
#   # Sentinel_shannon_poly <- as.polygons(Shannon, round=F, aggregate=F, extent=F,na.rm=F)
#   # writeVector(Sentinel_shannon_poly,filename=paste("~/data/output/INLA_modelling/Sentinel_shannon_poly", "res", window_size[x], "with_NA", sep="_"))
#   # writeVector(Arc_DEM_poly, filename = paste("~/data/output/INLA_modelling/Arc_DEM_poly","res",window_size[x], "with_NA", sep="_"))
# }
# 
# 
# ###############
# # Mean elevation
# ###############
# tile_DEM_proj <- rast("~/data/ArcDEM/tile_DEM_proj.tif")
# area_interest <- st_read("~/data/areas_of_interest.gpkg")
# area_interest_proj <- st_transform(area_interest,32613)
# ext(area_interest_proj)
# window_size <- c(10, 20, 30, 50, 100)
# multiple <- window_size*10
# extended_aoi <- ext(ext(area_interest_proj)[1]-multiple[5], ext(area_interest_proj)[2]+multiple[5], ext(area_interest_proj)[3]-multiple[5], ext(area_interest_proj)[4]+multiple[5])
# tile_DEM_crop <- crop(tile_DEM_proj, extended_aoi)
# fact <- window_size*10/2
# 
# for(x in 1:length(fact)){
#   # browser()
#   Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", window_size[x], "PC1278", sep="_"))
#   temp_rast <- rast(ext(Shannon), resolution = 2)
#   tile_DEM_crop_resample <- resample(tile_DEM_crop, temp_rast, method = "bilinear")
#   tile_DEM_crop_resample_noNA <- subst(tile_DEM_crop_resample, NA, -33) #have to adjust
#   tile_DEM_crop_aggregated <- aggregate(tile_DEM_crop_resample_noNA, fact=fact[x], fun="mean")
#   tile_DEM_masked <- mask(tile_DEM_crop_aggregated, Shannon)
#   writeRaster(tile_DEM_masked, filename = paste("~/data/ArcDEM/mean_ele_masked_", window_size[x], "_res.tif", sep=""))
#   # Arc_DEM_poly <- as.polygons(tile_DEM_masked, round=F, aggregate=F, extent=F, na.rm=F)
#   # Sentinel_shannon_poly <- as.polygons(Shannon, round=F, aggregate=F, extent=F,na.rm=F)
#   # writeVector(Sentinel_shannon_poly,filename=paste("~/data/output/INLA_modelling/Sentinel_shannon_poly", "res", window_size[x], "with_NA", sep="_"))
#   # writeVector(Arc_DEM_poly, filename = paste("~/data/output/INLA_modelling/Arc_DEM_poly","res",window_size[x], "with_NA", sep="_"))
# }
# 
# 
# ###############
# # Roughness
# ###############
# tile_DEM_proj <- rast("~/data/ArcDEM/tile_DEM_proj.tif")
# area_interest <- st_read("~/data/areas_of_interest.gpkg")
# area_interest_proj <- st_transform(area_interest,32613)
# ext(area_interest_proj)
# window_size <- c(10, 20, 30, 50, 100)
# multiple <- window_size*10
# extended_aoi <- ext(ext(area_interest_proj)[1]-multiple[5], ext(area_interest_proj)[2]+multiple[5], ext(area_interest_proj)[3]-multiple[5], ext(area_interest_proj)[4]+multiple[5])
# tile_DEM_crop <- crop(tile_DEM_proj, extended_aoi)
# fact <- window_size*10/2
# 
# for(x in 1:length(fact)){
#   # browser()
#   Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", window_size[x], "PC1278", sep="_"))
#   temp_rast <- rast(ext(Shannon), resolution = 2)
#   tile_DEM_crop_resample <- resample(tile_DEM_crop, temp_rast, method = "bilinear")
#   roughness <- terrain(tile_DEM_crop_resample, v="roughness")
#   roughness_noNA <- subst(roughness, NA, 0) #have to adjust
#   roughness_aggregated <- aggregate(roughness_noNA, fact=fact[x], fun="sd")
#   roughness_mask <- mask(roughness_aggregated, Shannon)
#   writeRaster(roughness_mask, filename = paste("~/data/ArcDEM/sd_roughness_masked_", window_size[x], "_res.tif", sep=""))
# }

