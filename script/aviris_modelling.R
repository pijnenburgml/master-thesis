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
  temp_rast <- rast(ext(Shannon), resolution=window_size[x]*5)  
  Shannon_resampled <- resample(Shannon, temp_rast)
  writeRaster(Shannon_resampled, filetype="ENVI", filename = paste("~/data/output/INLA_modelling/Aviris_modelling/Shannon_PC12_res_", window_size[x]*5, sep=""))
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

# sd(elevation)
for(x in 1:length(fact)){
  Shannon <- rast(paste("~/data/output/INLA_modelling/Aviris_modelling/Shannon_PC12_res_", window_size[x]*5, sep="")) ##have to add the name of the raster
  temp_rast <- rast(ext(Shannon), resolution = 2)
  tile_DEM_crop_resample <- resample(tile_DEM_crop, temp_rast, method = "bilinear")
  tile_DEM_crop_aggregated <- aggregate(tile_DEM_crop_resample, fact=fact[x], fun="sd")
  tile_DEM_masked <- mask(tile_DEM_crop_aggregated, Shannon)
  writeRaster(tile_DEM_masked, filename = paste("~/data/ArcDEM/sd_ele_aviris_", window_size[x], "_res.tif", sep=""))
}


# sd(slope)
for(x in 1:length(fact)){
  Shannon <- rast(paste("~/data/output/INLA_modelling/Aviris_modelling/Shannon_PC12_res_", window_size[x]*5, sep="")) ##have to add the name of the raster
  temp_rast <- rast(extended_aoi, resolution = 2)
  tile_DEM_crop_resample <- resample(tile_DEM_crop, temp_rast, method = "bilinear")
  slope <- terrain(tile_DEM_crop_resample, v="slope", neighbors=4, unit="degrees")
  slope_aggregated <- aggregate(slope, fact=fact[x], fun="sd")
  slope_aggregated_crop <- crop(slope_aggregated, Shannon)
  ext(slope_aggregated_crop) <- ext(Shannon)
  slope_mask <- mask(slope_aggregated_crop, Shannon)
  writeRaster(slope_mask, filename = paste("~/data/ArcDEM/sd_slope_aviris_", window_size[x], "_res.tif", sep=""), overwrite=T)
}


###########
# Modelling sd(slope)
###########
window_size <- c(20, 40, 60, 100, 200)
x <- 1
Shannon <- rast(paste("~/data/output/INLA_modelling/Aviris_modelling/Shannon_PC12_res_", window_size[x]*5, sep=""))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast(paste("~/data/ArcDEM/sd_slope_aviris_", window_size[x], "_res.tif", sep=""))
Slope_matrix <- as.matrix(Slope, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
sl <- inla.matrix2vector(Slope_matrix)
table(is.na(Slope_matrix))
table(is.na(Shannon_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, slope = sl, node=node)
data <- data[-c(which(is.na(s))),]

formula= Shannon_index ~ 1 + slope +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_slope_res_100 <- inla(formula,     
                                    family = "gamma",
                                    data = data,
                                    control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_slope_res_100)
save(Gamma_shannon_slope_res_100, file="~/data/output/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_100.Rdata")
get(load("~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_100.Rdata"))

# 200m
x <- 2
Shannon <- rast(paste("~/data/output/INLA_modelling/Aviris_modelling/Shannon_PC12_res_", window_size[x]*5, sep=""))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast(paste("~/data/ArcDEM/sd_slope_aviris_", window_size[x], "_res.tif", sep=""))
Slope_matrix <- as.matrix(Slope, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
sl <- inla.matrix2vector(Slope_matrix)
table(is.na(Slope_matrix))
table(is.na(Shannon_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, slope = sl, node=node)
data <- data[-c(which(is.na(s))),]

formula= Shannon_index ~ 1 + slope +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_slope_res_200 <- inla(formula,     
                                    family = "gamma",
                                    data = data,
                                    control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_slope_res_200)
save(Gamma_shannon_slope_res_200, file="~/data/output/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_200.Rdata")
get(load("~/scratch/INLA_modelling/Gamma_shannon_ele_slope_res_200.Rdata"))

# 300m
x <- 3
Shannon <- rast(paste("~/data/output/INLA_modelling/Aviris_modelling/Shannon_PC12_res_", window_size[x]*5, sep=""))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast(paste("~/data/ArcDEM/sd_slope_aviris_", window_size[x], "_res.tif", sep=""))
Slope_matrix <- as.matrix(Slope, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
sl <- inla.matrix2vector(Slope_matrix)
table(is.na(Slope_matrix))
table(is.na(Shannon_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, slope = sl, node=node)
data <- data[-c(which(is.na(s))),]

formula= Shannon_index ~ 1 + slope +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_slope_res_300 <- inla(formula,     
                                    family = "gamma",
                                    data = data,
                                    control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_slope_res_300)
save(Gamma_shannon_slope_res_300, file="~/data/output/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_300.Rdata")
get(load("~/scratch/INLA_modelling/Gamma_shannon_slope_res_300.Rdata"))


# 500
x <- 4
Shannon <- rast(paste("~/data/output/INLA_modelling/Aviris_modelling/Shannon_PC12_res_", window_size[x]*5, sep=""))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast(paste("~/data/ArcDEM/sd_slope_aviris_", window_size[x], "_res.tif", sep=""))
Slope_matrix <- as.matrix(Slope, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
sl <- inla.matrix2vector(Slope_matrix)
table(is.na(Slope_matrix))
table(is.na(Shannon_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, slope = sl, node=node)
data <- data[-c(which(is.na(s))),]

formula= Shannon_index ~ 1 + slope +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_slope_res_500 <- inla(formula,     
                                    family = "gamma",
                                    data = data,
                                    control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_slope_res_500)
save(Gamma_shannon_slope_res_500, file="~/data/output/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_500.Rdata")
get(load("~/scratch/INLA_modelling/Gamma_shannon_slope_res_500.Rdata"))


# 1000m
x <- 5
Shannon <- rast(paste("~/data/output/INLA_modelling/Aviris_modelling/Shannon_PC12_res_", window_size[x]*5, sep=""))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast(paste("~/data/ArcDEM/sd_slope_aviris_", window_size[x], "_res.tif", sep=""))
Slope_matrix <- as.matrix(Slope, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
sl <- inla.matrix2vector(Slope_matrix)
table(is.na(Slope_matrix))
table(is.na(Shannon_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, slope = sl, node=node)
data <- data[-c(which(is.na(s))),]

formula= Shannon_index ~ 1 + slope +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_slope_res_1000 <- inla(formula,     
                                     family = "gamma",
                                     data = data,
                                     control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_slope_res_1000)
save(Gamma_shannon_slope_res_1000, file="~/data/output/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_1000.Rdata")
get(load("~/scratch/INLA_modelling/Gamma_shannon_slope_res_1000.Rdata"))



###########
# Modelling sd(ele)
###########
window_size <- c(20, 40, 60, 100, 200)
x <- 1
Shannon <- rast(paste("~/data/output/INLA_modelling/Aviris_modelling/Shannon_PC12_res_", window_size[x]*5, sep=""))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Ele <- rast(paste("~/data/ArcDEM/sd_ele_aviris_", window_size[x], "_res.tif", sep=""))
Ele_matrix <- as.matrix(Ele, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
e <- inla.matrix2vector(Ele_matrix)
table(is.na(Ele_matrix))
table(is.na(Shannon_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, Elevation = e, node=node)
data <- data[-c(which(is.na(s))),]

formula= Shannon_index ~ 1 + Elevation +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_ele_res_100 <- inla(formula,     
                                    family = "gamma",
                                    data = data,
                                    control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_ele_res_100)
save(Gamma_shannon_ele_res_100, file="~/data/output/INLA_modelling/Aviris_modelling/Gamma_shannon_ele_res_100.Rdata")
get(load("~/data/output/INLA_modelling/Aviris_modelling/Gamma_shannon_ele_res_100.Rdata"))

# 200m 
x <- 2
Shannon <- rast(paste("~/data/output/INLA_modelling/Aviris_modelling/Shannon_PC12_res_", window_size[x]*5, sep=""))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Ele <- rast(paste("~/data/ArcDEM/sd_ele_aviris_", window_size[x], "_res.tif", sep=""))
Ele_matrix <- as.matrix(Ele, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
e <- inla.matrix2vector(Ele_matrix)
table(is.na(Ele_matrix))
table(is.na(Shannon_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, Elevation = e, node=node)
data <- data[-c(which(is.na(s))),]

formula= Shannon_index ~ 1 + Elevation +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_ele_res_200 <- inla(formula,     
                                  family = "gamma",
                                  data = data,
                                  control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_ele_res_200)
save(Gamma_shannon_ele_res_200, file="~/data/output/INLA_modelling/Aviris_modelling/Gamma_shannon_ele_res_200.Rdata")
get(load("~/data/output/INLA_modelling/Aviris_modelling/Gamma_shannon_ele_res_200.Rdata"))

# 300m 
x <- 3
Shannon <- rast(paste("~/data/output/INLA_modelling/Aviris_modelling/Shannon_PC12_res_", window_size[x]*5, sep=""))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Ele <- rast(paste("~/data/ArcDEM/sd_ele_aviris_", window_size[x], "_res.tif", sep=""))
Ele_matrix <- as.matrix(Ele, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
e <- inla.matrix2vector(Ele_matrix)
table(is.na(Ele_matrix))
table(is.na(Shannon_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, Elevation = e, node=node)
data <- data[-c(which(is.na(s))),]

formula= Shannon_index ~ 1 + Elevation +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_ele_res_300 <- inla(formula,     
                                  family = "gamma",
                                  data = data,
                                  control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_ele_res_300)
save(Gamma_shannon_ele_res_300, file="~/data/output/INLA_modelling/Aviris_modelling/Gamma_shannon_ele_res_300.Rdata")
get(load("~/data/output/INLA_modelling/Aviris_modelling/Gamma_shannon_ele_res_300.Rdata"))


# 500m 
x <- 4
Shannon <- rast(paste("~/data/output/INLA_modelling/Aviris_modelling/Shannon_PC12_res_", window_size[x]*5, sep=""))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Ele <- rast(paste("~/data/ArcDEM/sd_ele_aviris_", window_size[x], "_res.tif", sep=""))
Ele_matrix <- as.matrix(Ele, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
e <- inla.matrix2vector(Ele_matrix)
table(is.na(Ele_matrix))
table(is.na(Shannon_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, Elevation = e, node=node)
data <- data[-c(which(is.na(s))),]

formula= Shannon_index ~ 1 + Elevation +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_ele_res_500 <- inla(formula,     
                                  family = "gamma",
                                  data = data,
                                  control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_ele_res_500)
save(Gamma_shannon_ele_res_500, file="~/data/output/INLA_modelling/Aviris_modelling/Gamma_shannon_ele_res_500.Rdata")
get(load("~/data/output/INLA_modelling/Aviris_modelling/Gamma_shannon_ele_res_500.Rdata"))


# 1000m 
x <- 5
Shannon <- rast(paste("~/data/output/INLA_modelling/Aviris_modelling/Shannon_PC12_res_", window_size[x]*5, sep=""))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Ele <- rast(paste("~/data/ArcDEM/sd_ele_aviris_", window_size[x], "_res.tif", sep=""))
Ele_matrix <- as.matrix(Ele, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
e <- inla.matrix2vector(Ele_matrix)
table(is.na(Ele_matrix))
table(is.na(Shannon_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, Elevation = e, node=node)
data <- data[-c(which(is.na(s))),]

formula= Shannon_index ~ 1 + Elevation +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_ele_res_1000 <- inla(formula,     
                                  family = "gamma",
                                  data = data,
                                  control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_ele_res_1000)
save(Gamma_shannon_ele_res_1000, file="~/data/output/INLA_modelling/Aviris_modelling/Gamma_shannon_ele_res_1000.Rdata")
get(load("~/data/output/INLA_modelling/Aviris_modelling/Gamma_shannon_ele_res_1000.Rdata"))












############
# prepare corresponding Sentinel data
############
Datadir <- "~/data/biodivmapR_sent"
NameRaster <- "sent_crop_envi_BIL"
Input_Image_File <- file.path(Datadir,NameRaster)
NameMask <- "mask_matchin_aviris"
Input_Mask_File <- file.path(Datadir, NameMask)
Output_Dir <- "~/data/biodivmapR_sent/RESULTS_aviris_extend"
Continuum_Removal <- F
TypePCA <- 'SPCA'
FilterPCA <- F
nbCPU <- 4
MaxRAM <- 4
nbclusters <- 20
#PCA output
PCA_Output <- get(load("~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/PCA/PCA_Output.Rdata"))
SelectedPCs = c(1,2,7,8)
#Clustering output
Kmeans_info <- get(load("~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/SpectralSpecies/Kmeans_info_PC1278.Rdata")) #have to adjust name
window_size <- c(10, 20, 30, 50, 100)
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


#############
# preparing elevation data
#############
tile_DEM_proj <- rast("~/data/ArcDEM/tile_DEM_proj.tif")
area_interest <- st_read("~/data/areas_of_interest.gpkg")
area_interest_proj <- st_transform(area_interest,32613)
ext(area_interest_proj)
window_size <- c(10, 20, 30, 50, 100)
multiple <- window_size*10
extended_aoi <- ext(ext(area_interest_proj)[1]-multiple[5], ext(area_interest_proj)[2]+multiple[5], ext(area_interest_proj)[3]-multiple[5], ext(area_interest_proj)[4]+multiple[5])
tile_DEM_crop <- crop(tile_DEM_proj, extended_aoi)
fact <- window_size*10/2

# sd(slope)
for(x in 4:length(fact)){
  Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_", window_size[x], "_PC1278", sep="")) ##have to add the name of the raster
  temp_rast <- rast(extended_aoi, resolution = 2)
  tile_DEM_crop_resample <- resample(tile_DEM_crop, temp_rast, method = "bilinear")
  slope <- terrain(tile_DEM_crop_resample, v="slope", neighbors=4, unit="degrees")
  slope_aggregated <- aggregate(slope, fact=fact[x], fun="sd")
  slope_aggregated_crop <- crop(slope_aggregated, Shannon)
  ext(slope_aggregated_crop) <- ext(Shannon)
  slope_mask <- mask(slope_aggregated_crop, Shannon)
  writeRaster(slope_mask, filename = paste("~/data/ArcDEM/sd_slope_sent_aviris_aoi_", window_size[x], "_res.tif", sep=""), overwrite=T)
}

tile_DEM_crop <- crop(tile_DEM_proj, area_interest_proj)
# sd(elevation)
for(x in 1:length(fact)){
  Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_", window_size[x], "_PC1278", sep="")) ##have to add the name of the raster
  temp_rast <- rast(extended_aoi, resolution = 2)
  tile_DEM_crop_resample <- resample(tile_DEM_crop, temp_rast, method = "bilinear")
  tile_DEM_crop_aggregated <- aggregate(tile_DEM_crop_resample, fact=fact[x], fun="sd")
  tile_DEM_crop_crop <- crop(tile_DEM_crop_aggregated, Shannon)
  ext(tile_DEM_crop_crop) <- ext(Shannon)
  tile_DEM_masked <- mask(tile_DEM_crop_crop, Shannon)
  writeRaster(tile_DEM_masked, filename = paste("~/data/ArcDEM/sd_ele_sent_aviris_aoi_", window_size[x], "_res.tif", sep=""))
}



###########
# Modelling sd(slope)
###########

# 100m 
window_size <- c(10, 20, 30, 50, 100)
x <- 1
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_", window_size[x], "_PC1278", sep=""))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast(paste("~/data/ArcDEM/sd_slope_sent_aviris_aoi_", window_size[x], "_res.tif", sep=""))
Slope_matrix <- as.matrix(Slope, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
sl <- inla.matrix2vector(Slope_matrix)
table(is.na(Slope_matrix))
table(is.na(Shannon_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, slope = sl, node=node)
data <- data[-c(which(is.na(s))),]

formula= Shannon_index ~ 1 + slope +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_slope_res_100 <- inla(formula,     
                                    family = "gamma",
                                    data = data,
                                    control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_slope_res_100)
save(Gamma_shannon_slope_res_100, file="~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_100.Rdata")
get(load("~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_100.Rdata"))

# 200m 
window_size <- c(10, 20, 30, 50, 100)
x <- 2
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_", window_size[x], "_PC1278", sep=""))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast(paste("~/data/ArcDEM/sd_slope_sent_aviris_aoi_", window_size[x], "_res.tif", sep=""))
Slope_matrix <- as.matrix(Slope, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
sl <- inla.matrix2vector(Slope_matrix)
table(is.na(Slope_matrix))
table(is.na(Shannon_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, slope = sl, node=node)
data <- data[-c(which(is.na(s))),]

formula= Shannon_index ~ 1 + slope +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_slope_res_200 <- inla(formula,     
                                    family = "gamma",
                                    data = data,
                                    control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_slope_res_200)
save(Gamma_shannon_slope_res_200, file="~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_200.Rdata")
get(load("~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_200.Rdata"))


# 300m 
window_size <- c(10, 20, 30, 50, 100)
x <- 3
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_", window_size[x], "_PC1278", sep=""))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast(paste("~/data/ArcDEM/sd_slope_sent_aviris_aoi_", window_size[x], "_res.tif", sep=""))
Slope_matrix <- as.matrix(Slope, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
sl <- inla.matrix2vector(Slope_matrix)
table(is.na(Slope_matrix))
table(is.na(Shannon_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, slope = sl, node=node)
data <- data[-c(which(is.na(s))),]

formula= Shannon_index ~ 1 + slope +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_slope_res_300 <- inla(formula,     
                                    family = "gamma",
                                    data = data,
                                    control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_slope_res_300)
save(Gamma_shannon_slope_res_300, file="~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_300.Rdata")
get(load("~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_300.Rdata"))


# 500m 
window_size <- c(10, 20, 30, 50, 100)
x <- 4
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_", window_size[x], "_PC1278", sep=""))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast(paste("~/data/ArcDEM/sd_slope_sent_aviris_aoi_", window_size[x], "_res.tif", sep=""))
Slope_matrix <- as.matrix(Slope, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
sl <- inla.matrix2vector(Slope_matrix)
table(is.na(Slope_matrix))
table(is.na(Shannon_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, slope = sl, node=node)
data <- data[-c(which(is.na(s))),]

formula= Shannon_index ~ 1 + slope +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_slope_res_500 <- inla(formula,     
                                    family = "gamma",
                                    data = data,
                                    control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_slope_res_500)
save(Gamma_shannon_slope_res_500, file="~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_500.Rdata")
get(load("~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_500.Rdata"))

# 1000m 
window_size <- c(10, 20, 30, 50, 100)
x <- 5
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_", window_size[x], "_PC1278", sep=""))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast(paste("~/data/ArcDEM/sd_slope_sent_aviris_aoi_", window_size[x], "_res.tif", sep=""))
Slope_matrix <- as.matrix(Slope, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
sl <- inla.matrix2vector(Slope_matrix)
table(is.na(Slope_matrix))
table(is.na(Shannon_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, slope = sl, node=node)
data <- data[-c(which(is.na(s))),]

formula= Shannon_index ~ 1 + slope +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_slope_res_1000 <- inla(formula,     
                                    family = "gamma",
                                    data = data,
                                    control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_slope_res_1000)
save(Gamma_shannon_slope_res_1000, file="~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_1000.Rdata")
get(load("~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_1000.Rdata"))



###########
# Modelling sd(elevation)
###########

# 100m 
window_size <- c(10, 20, 30, 50, 100)
x <- 1
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_", window_size[x], "_PC1278", sep=""))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Elevation <- rast(paste("~/data/ArcDEM/sd_slope_sent_aviris_aoi_", window_size[x], "_res.tif", sep=""))
Ele_matrix <- as.matrix(Elevation, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
e <- inla.matrix2vector(Ele_matrix)
table(is.na(Ele_matrix))
table(is.na(Shannon_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, Ele = e, node=node)
data <- data[-c(which(is.na(s))),]

formula= Shannon_index ~ 1 + Ele +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_ele_res_100 <- inla(formula,     
                                    family = "gamma",
                                    data = data,
                                    control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_ele_res_100)
save(Gamma_shannon_ele_res_100, file="~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_ele_res_100.Rdata")
get(load("~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_ele_res_100.Rdata"))

# 200m 
window_size <- c(10, 20, 30, 50, 100)
x <- 2
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_", window_size[x], "_PC1278", sep=""))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Elevation <- rast(paste("~/data/ArcDEM/sd_slope_sent_aviris_aoi_", window_size[x], "_res.tif", sep=""))
Ele_matrix <- as.matrix(Elevation, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
sl <- inla.matrix2vector(Ele_matrix)
table(is.na(Ele_matrix))
table(is.na(Shannon_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, slope = sl, node=node)
data <- data[-c(which(is.na(s))),]

formula= Shannon_index ~ 1 + slope +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_ele_res_200 <- inla(formula,     
                                  family = "gamma",
                                  data = data,
                                  control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_ele_res_200)
save(Gamma_shannon_ele_res_200, file="~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_ele_res_200.Rdata")
get(load("~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_ele_res_200.Rdata"))


# 300m 
window_size <- c(10, 20, 30, 50, 100)
x <- 3
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_", window_size[x], "_PC1278", sep=""))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Elevation <- rast(paste("~/data/ArcDEM/sd_slope_sent_aviris_aoi_", window_size[x], "_res.tif", sep=""))
Ele_matrix <- as.matrix(Elevation, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
sl <- inla.matrix2vector(Ele_matrix)
table(is.na(Ele_matrix))
table(is.na(Shannon_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, slope = sl, node=node)
data <- data[-c(which(is.na(s))),]

formula= Shannon_index ~ 1 + slope +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_ele_res_300 <- inla(formula,     
                                  family = "gamma",
                                  data = data,
                                  control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_ele_res_300)
save(Gamma_shannon_ele_res_300, file="~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_ele_res_300.Rdata")
get(load("~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_ele_res_300.Rdata"))


# 500m 
window_size <- c(10, 20, 30, 50, 100)
x <- 4
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_", window_size[x], "_PC1278", sep=""))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Elevation <- rast(paste("~/data/ArcDEM/sd_slope_sent_aviris_aoi_", window_size[x], "_res.tif", sep=""))
Ele_matrix <- as.matrix(Elevation, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
sl <- inla.matrix2vector(Ele_matrix)
table(is.na(Ele_matrix))
table(is.na(Shannon_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, slope = sl, node=node)
data <- data[-c(which(is.na(s))),]

formula= Shannon_index ~ 1 + slope +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_ele_res_500 <- inla(formula,     
                                  family = "gamma",
                                  data = data,
                                  control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_ele_res_500)
save(Gamma_shannon_ele_res_500, file="~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_ele_res_500.Rdata")
get(load("~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_ele_res_500.Rdata"))

# 1000m 
window_size <- c(10, 20, 30, 50, 100)
x <- 5
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_", window_size[x], "_PC1278", sep=""))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Elevation <- rast(paste("~/data/ArcDEM/sd_slope_sent_aviris_aoi_", window_size[x], "_res.tif", sep=""))
Ele_matrix <- as.matrix(Elevation, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
sl <- inla.matrix2vector(Ele_matrix)
table(is.na(Ele_matrix))
table(is.na(Shannon_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, slope = sl, node=node)
data <- data[-c(which(is.na(s))),]

formula= Shannon_index ~ 1 + slope +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_ele_res_100 <- inla(formula,     
                                  family = "gamma",
                                  data = data,
                                  control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_ele_res_1000)
save(Gamma_shannon_ele_res_1000, file="~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_ele_res_1000.Rdata")
get(load("~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_ele_res_1000.Rdata"))








