#########
# Modelling with aviris data
#########
setwd("~/data")
library(sf)
library(tidyverse)
library(dplyr)
library(tidyterra)
library(remotes)
library(dissUtils)
library(biodivMapR)
library(ggplot2)
library(cowplot)

library(sp)
library(rgdal)  
require(spdep)
library(INLA)
library(terra)
library(INLAutils)
library(scales)



#############
# preparing Shannon index data
#############
Datadir <- "/scratch/mpijne/reflectance_data"
NameRaster <- "ang20190802t220708_rfl_rect"
Input_Image_File <- file.path(Datadir,NameRaster)
NameMask <- "strip_0708_aoi_mask"
Input_Mask_File <- file.path(Datadir, NameMask)
Output_Dir <- "biodivmapR_Aviris/RESULTS"
Continuum_Removal <- F
TypePCA <- 'SPCA'
FilterPCA <- F
window_size <- 20
nbCPU <- 8
MaxRAM <- 10
nbclusters <- 20

#PCA output
PCA_Output <- get(load("~/data/biodivmapR_Aviris/RESULTS/ang20190802t220708_rfl_rect/SPCA/PCA/PCA_Output.Rdata"))
SelectedPCs = c(1,2)

#Clustering output
Kmeans_info <- get(load("~/data/biodivmapR_Aviris/RESULTS/ang20190802t220708_rfl_rect/SPCA/SpectralSpecies/Kmeans_info_PC12.Rdata")) #have to adjust name
window_size <- c(20, 40, 60, 100, 200)
Index_Alpha <- c("Shannon")

for (x in 1:length(window_size)){	
  map_alpha_div_ML(Input_Image_File = Input_Image_File,
                   Output_Dir = Output_Dir,
                   TypePCA = TypePCA,
                   window_size = window_size[x],
                   nbCPU = nbCPU,
                   MaxRAM = MaxRAM,
                   Index_Alpha = Index_Alpha,
                   nbclusters = nbclusters, SelectedPCs = SelectedPCs)
}

# result in maps of resolution 100m, 200m, 300m, 500m, 1000m.

for (x in 1:length(window_size)) {
  Shannon <- rast(paste("~/data/biodivmapR_Aviris/RESULTS/ang20190802t220708_rfl_rect/SPCA/ALPHA/Shannon_", window_size[x], "_PC12", sep=""))
  temp_rast <- rast(ext(Shannon, resolution=window_size[x]*5))  
  Shannon_resampled <- resample(Shannon, temp_rast)
  writeRaster(Shannon_resampled, filename = paste("~/data/output/INLA_modelling/Aviris_modelling/Shannon_PC12_res_", window_size[x]*5, sep=""))
}


#############
# preparing elevation data
#############
tile_DEM_proj <- rast("~/data/ArcDEM/tile_DEM_proj.tif")
cell <- "ang20190802t220708_rfl_rect"
tile <- rast(file.path("/scratch/mpijne/reflectance_data", cell))
ext(tile)
window_size <- c(20, 40, 60, 100, 200)
fact <- c(50, 100, 150, 250, 500)
multiple <- window_size*10
extended_aoi <- ext(ext(tile)[1]-multiple[5], ext(tile)[2]+multiple[5], ext(tile)[3]-multiple[5], ext(tile)[4]+multiple[5])
tile_DEM_crop <- crop(tile_DEM_proj, extended_aoi)

for(x in 1:length(fact)){
  # browser()
  Shannon <- rast("~/data/output/INLA_modelling/Aviris_modelling/") ##have to add the name of the raster
  temp_rast <- rast(ext(Shannon), resolution = 2)
  tile_DEM_crop_resample <- resample(tile_DEM_crop, temp_rast, method = "bilinear")
  tile_DEM_crop_aggregated <- aggregate(tile_DEM_crop_resample, fact=fact[x], fun="sd")
  tile_DEM_masked <- mask(tile_DEM_crop_aggregated, Shannon)
  writeRaster(tile_DEM_masked, filename = paste("~/data/ArcDEM/sd_ele_agg10_", window_size[x], "_res.tif", sep=""))
}






















