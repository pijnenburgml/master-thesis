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


#################
# Scaling analysis 
#################

################################################################################
Datadir <- "~/data/biodivmapR_sent"
NameRaster <- "sent_crop_envi_BIL"
# Define path for image file to be processed
Input_Image_File <- file.path(Datadir,NameRaster)
# Define path for corresponding mask file
NameMask <- "mask_sent2_final_NA"
Input_Mask_File <- file.path(Datadir, NameMask)
# Input_Mask_File <- F
# Define path for master output directory where files produced during the process are saved
Output_Dir <- '~/data/biodivmapR_sent/RESULTS_cluster_20'
# dir.create(path = Output_Dir,recursive = T,showWarnings = F)
# NDVI_Thresh <- 0.8
# Blue_Thresh <- 500
# NIR_Thresh <- 1500
# Apply normalization with continuum removal?
Continuum_Removal <- F
# Type of dimensionality reduction
TypePCA <- 'SPCA'
# PCA FILTERING:        Set to TRUE if you want second filtering based on PCA outliers to be processed.
# Slower process
# Automatically set to FALSE if TypePCA     = 'MNF'
FilterPCA <- F
# window size forcomputation of spectral diversity
# window_size <-10
# # computational parameters
nbCPU <- 2
MaxRAM <- 12
# number of clusters (spectral species)
nbclusters <- 20


################################################################################
##                Perform alpha and beta diversity mapping                    ##
## https://jbferet.github.io/biodivMapR/articles/biodivMapR_6.html            ##
################################################################################

# To have the name of the PCs selected to map the alpha diversity in the filename 
# run the function saved in the document map_alpha_div_PC_naming.R to have it in the 
# global environment

#PCA output
PCA_Output <- get(load("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/PCA/PCA_Output.RData"))
SelectedPCs = c(1,2,7,8)

#Clustering output
Kmeans_info <- get(load("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/SpectralSpecies/Kmeans_info_PC1278.Rdata")) #have to adjust name
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

# result in maps of resolution 100m, 200m, 300m, 500m, 1000m. 

#######
# assessment
#######
viridis_colors <- viridis::plasma(20)
for(x in 1:length(window_size)){
  path <- paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_", window_size[x],"_PC1278", sep="")
  m <- rast(path)
  plot(m, col=viridis_colors)
    
}

#####
# prepare elevation data
#####
tile_DEM_proj <- rast("~/data/ArcDEM/tile_DEM_proj.tif")
area_interest <- st_read("~/data/areas_of_interest.gpkg")
area_interest_proj <- st_transform(area_interest,32613)
ext(area_interest_proj)
multiple <- window_size*10
extended_aoi <- ext(ext(area_interest_proj)[1]-multiple[5], ext(area_interest_proj)[2]+multiple[5], ext(area_interest_proj)[3]-multiple[5], ext(area_interest_proj)[4]+multiple[5])
tile_DEM_crop <- crop(tile_DEM_proj, extended_aoi)

window_size <- c(10, 20, 30, 50, 100)
fact <- window_size*10/2

for(x in 1:length(fact)){
  # browser()
  Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", window_size[x], "PC1278", sep="_"))
  temp_rast <- rast(ext(Shannon), resolution = 2)
  tile_DEM_crop_resample <- resample(tile_DEM_crop, temp_rast, method = "bilinear")
  tile_DEM_crop_resample_noNA <- subst(tile_DEM_crop_resample, NA, -33) #have to adjust
  tile_DEM_crop_aggregated <- aggregate(tile_DEM_crop_resample_noNA, fact=fact[x], fun="sd")
  tile_DEM_masked <- mask(tile_DEM_crop_aggregated, Shannon)
  writeRaster(tile_DEM_masked, filename = paste("~/data/ArcDEM/ArcDEM_masked_", window_size[x], "_res.tif", sep=""))
  # Arc_DEM_poly <- as.polygons(tile_DEM_masked, round=F, aggregate=F, extent=F, na.rm=F)
  # Sentinel_shannon_poly <- as.polygons(Shannon, round=F, aggregate=F, extent=F,na.rm=F)
  # writeVector(Sentinel_shannon_poly,filename=paste("~/data/output/INLA_modelling/Sentinel_shannon_poly", "res", window_size[x], "with_NA", sep="_"))
  # writeVector(Arc_DEM_poly, filename = paste("~/data/output/INLA_modelling/Arc_DEM_poly","res",window_size[x], "with_NA", sep="_"))
}

# x <- 1 #to be change!
# path_sent <- paste("~/data/output/INLA_modelling/Sentinel_shannon_poly_res_", window_size[x],"_with_NA", "/Sentinel_shannon_poly_res_", window_size[x], "_with_NA.shp", sep="")
# path_arcdem <- paste("~/data/output/INLA_modelling/Arc_DEM_poly_res_", window_size[x], "_with_NA", "/Arc_DEM_poly_res_", window_size[x], "_with_NA.shp", sep="")
# Sentinel_shannon_vect <- vect(path_sent)
# Arc_DEM_vect <- vect(path_arcdem)
# sd_topo <- Arc_DEM_vect$X29_21_1_1
# Sentinel_shannon_vect$sd_topo <- sd_topo
# writeVector(Sentinel_shannon_vect, filename = paste("~/data/output/INLA_modelling/model_object_res_", window_size[x], "_with_NA.shp", sep=""), overwrite=T)

#################
# Modelling
#################

Sentinel_lattice <-readOGR("~/data/output/INLA_modelling/model_object_res_10_with_NA.shp")
Sentinel_data <- Sentinel_lattice@data
Sentinel_data$Shannon_10[!is.na(Sentinel_data$Shannon_10)] <- Sentinel_data$Shannon_10[!is.na(Sentinel_data$Shannon_10)]+ 1
which(Sentinel_data$Shannon_10==0)
hist(Sentinel_data$Shannon_10)
hist(Sentinel_data$sd_topo)
hist(log(Sentinel_data$sd_topo))
plot(Sentinel_data$Shannon_10~Sentinel_data$sd_topo)
plot(Sentinel_data$Shannon_10~log(Sentinel_data$sd_topo))

Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 10, "PC1278", sep="_"))
Elev <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 10, "_res.tif", sep=""))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Elev_matrix <- as.matrix(Elev, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
e <- inla.matrix2vector(Elev_matrix)
table(is.na(Elev_matrix))
node = 1:n
node[which(is.na(s))] <- NA
table(is.na(node))
data <- data.frame(Shannon_index = s, topo = e, spatial=node)
y_noNA <- y[-c(which(is.na(y)))]

# node = 1:length(y_noNA)
# Sentinel_data$node <- node

# log.range = list(initial = log(5), fixed=TRUE) 
# hyperpar_matern = list(initial = -3, param=c(23.36,0.001))
# formula_matern = y ~ wildfire + logging + lichwood + openlich + deciduous + 
#   water + wetland + meanelev +
#   f(node_matern, model = "matern2d", nrow = nrow.larger, 
#     ncol = ncol.larger, hyper = list(range = log.range, prec 
#                                      = hyperpar_matern))
# hyper = list(range = list(param =c(1, 1),prior = "loggamma",initial=1),
#              prec = list(param=c(1, 1)))

formula= Shannon_10 ~ 1+ log(sd_topo)+
  f(node, model="matern2d", nrow=nrow, ncol=ncol,
    hyper = list(range = list(initual=log(1000))))

model_lattice_gaussian_loglink_withNA <- inla(formula,     
                                       family = "gaussian",
                                       control.family=list(link='log'),
                                       data = Sentinel_data,
                                       control.compute = list(cpo = T, dic = T, waic = T, return.marginals.predictor=TRUE), verbose=TRUE)


summary(model_lattice_gaussian_loglink)


nrow=20
ncol=30
n = nrow*ncol
s.noise = 1
zi.mat = matrix(NA,nrow=nrow,ncol=ncol)
i=1:nrow
for(j in 1:ncol)
  zi.mat[i,j] = 3*exp(-(i-j)^2/4)
## iid noise
noise.mat=matrix(rnorm(nrow*ncol, sd=s.noise),nrow,ncol)
## make simulated data with no spatial component
y.mat = zi.mat + noise.mat
plot(rast(y.mat))
y.mat[15,4] <- NA
y.mat[3,24] <- NA
plot(rast(y.mat))
## convert matrices to the internal representation in INLA
y = inla.matrix2vector(y.mat)
node = 1:n
formula= y ~ 1+ f(node, model="matern2d", nu=1, nrow=nrow, ncol=ncol,
                  hyper = list(range = list(param =c(1, 1),
                                            prior = "loggamma",initial=1),
                               prec = list(param=c(1, 1))))

data=data.frame(y=y,node=node)
result=inla(formula, family="gaussian", data=data, verbose=TRUE,
            control.predictor = list(compute = TRUE),
            control.family = list(hyper = list(theta = list(initial = log(1/s.noise^2),
                                                            fixed = FALSE))),
            keep=T)



