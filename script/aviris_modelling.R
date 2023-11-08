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
library(grid)
library(gridExtra)

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
nbCPU <- 20
MaxRAM <- 8
nbclusters <- 20

#PCA output
PCA_Output <- get(load("~/data/biodivmapR_Aviris/RESULTS/ang20190802t220708_rfl_rect/SPCA/PCA/PCA_Output.Rdata"))
SelectedPCs = c(1,2,8,9)

#Clustering output
Kmeans_info <- get(load("~/data/biodivmapR_Aviris/RESULTS/ang20190802t220708_rfl_rect/SPCA/SpectralSpecies/Kmeans_info_PC1289.Rdata")) #have to adjust name
window_size <- c(20, 40, 60)
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

# result in maps of resolution 100m, 200m, 300m.


#############
# preparing elevation data
#############
tile_DEM_proj <- rast("~/data/ArcDEM/tile_DEM_proj.tif")
cell <- "ang20190802t220708_rfl_rect"
tile <- rast(file.path("/scratch/mpijne/reflectance_data", cell))
ext(tile)
window_size <- c(20, 40, 60)
fact <- c(50, 100, 150)
multiple <- window_size*10
extended_aoi <- ext(ext(tile)[1]-multiple[3], ext(tile)[2]+multiple[3], ext(tile)[3]-multiple[3], ext(tile)[4]+multiple[3])
tile_DEM_crop <- crop(tile_DEM_proj, extended_aoi)


# sd(slope)
for(x in 1:length(fact)){
  Shannon <- rast(paste("~/data/biodivmapR_Aviris/RESULTS/ang20190802t220708_rfl_rect/SPCA/ALPHA/Shannon_", window_size[x], "_PC1289", sep="")) ##have to add the name of the raster
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
window_size <- c(20, 40, 60)
x <- 1
Shannon <- rast(paste("~/data/biodivmapR_Aviris/RESULTS/ang20190802t220708_rfl_rect/SPCA/ALPHA/Shannon_", window_size[x], "_PC1289", sep="")) ##have to add the name of the raster
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast(paste("~/data/ArcDEM/sd_slope_aviris_", window_size[x], "_res.tif", sep=""))
Slope_matrix <- as.matrix(Slope, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
length(which(s==0))
s[!is.na(s)] <- s[!is.na(s)]+ 1
sl <- inla.matrix2vector(Slope_matrix)
table(is.na(Slope_matrix))
table(is.na(Shannon_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, slope = sl, node=node)
data <- data[-c(which(is.na(s))),]
nrow(data)

formula= Shannon_index ~ 1 + slope +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_slope_res_100 <- inla(formula,
                                    family = "gamma",
                                    data = data,
                                    control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_slope_res_100)
save(Gamma_shannon_slope_res_100, file="~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_100.Rdata")
get(load("~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_100.Rdata"))

observed <- data$Shannon_index
fit <- Gamma_shannon_slope_res_100$summary.fitted.values$mean[1:length(observed)]
pseudo_r2=function(observed,fit){
  res =  observed-fit
  RRes=sum((res)^2,na.rm = T)
  RRtot=sum((observed-mean(fit,na.rm=T))^2,na.rm = T)
  pseudo_r2_val=1-RRes/RRtot
  print(RRes)
  print(RRtot)
  return(pseudo_r2_val)  
}
pr <- round(pseudo_r2(observed, fit), digits = 3)
qq100 <- ggplot()+  
  geom_abline(slope=1, intercept=0, col="grey")+
  geom_point(data=data.frame(observed=observed, fit=fit), aes(fit, observed), pch=21)+
  ylab("")+
  xlab("")+
  scale_x_continuous(breaks=c(1,2,3), limits =c(1, 3.5))+
  scale_y_continuous(breaks=c(1,2,3), limits =c(1, 4))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
        panel.grid.major.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=margin(t = 1, r = 100, b = 1, l = 1, unit = "pt"), 
        panel.background = element_rect(fill=NA),
        text = element_text(size=20))+
  draw_text(text=paste("pseudo-R2 =", pr), x=1.7, y=3, size=20)
qq100
-sum(log(Gamma_shannon_slope_res_100$cpo$cpo))


# 200m
x <- 2
Shannon <- rast(paste("~/data/biodivmapR_Aviris/RESULTS/ang20190802t220708_rfl_rect/SPCA/ALPHA/Shannon_", window_size[x], "_PC1289", sep="")) ##have to add the name of the raster
zero_color <- "red"
plot(Shannon$Shannon_40_PC1289, col = c(zero_color, grey.colors(30)))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast(paste("~/data/ArcDEM/sd_slope_aviris_", window_size[x], "_res.tif", sep=""))
Slope_matrix <- as.matrix(Slope, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
length(which(s==0))
s[!is.na(s)] <- s[!is.na(s)]+ 1
sl <- inla.matrix2vector(Slope_matrix)
table(is.na(Slope_matrix))
table(is.na(Shannon_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, slope = sl, node=node)
data <- data[-c(which(is.na(s))),]
nrow(data)

formula= Shannon_index ~ 1 + slope +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_slope_res_200 <- inla(formula,
                                    family = "gamma",
                                    data = data,
                                    control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_slope_res_200)
save(Gamma_shannon_slope_res_200, file="~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_200.Rdata")
get(load("~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_200.Rdata"))

observed <- data$Shannon_index
fit <- Gamma_shannon_slope_res_200$summary.fitted.values$mean[1:length(observed)]
pseudo_r2=function(observed,fit){
  res =  observed-fit
  RRes=sum((res)^2,na.rm = T)
  RRtot=sum((observed-mean(fit,na.rm=T))^2,na.rm = T)
  pseudo_r2_val=1-RRes/RRtot
  print(RRes)
  print(RRtot)
  return(pseudo_r2_val)  
}
pr <- round(pseudo_r2(observed, fit), digits = 3)
qq200 <- ggplot()+  
  geom_abline(slope=1, intercept=0, col="grey")+
  geom_point(data=data.frame(observed=observed, fit=fit), aes(fit, observed), pch=21)+
  ylab("")+
  xlab("")+
  scale_x_continuous(breaks=c(1,2,3), limits =c(1, 3.5))+
  scale_y_continuous(breaks=c(1,2,3), limits =c(1, 4))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
        panel.grid.major.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=margin(t = 1, r = 100, b = 1, l = 1, unit = "pt"), 
        panel.background = element_rect(fill=NA),
        text = element_text(size=20))+
  draw_text(text=paste("pseudo-R2 =", pr), x=1.7, y=3, size=20)
qq200
-sum(log(Gamma_shannon_slope_res_200$cpo$cpo))


# 300m
x <- 3
Shannon <- rast(paste("~/data/biodivmapR_Aviris/RESULTS/ang20190802t220708_rfl_rect/SPCA/ALPHA/Shannon_", window_size[x], "_PC1289", sep="")) ##have to add the name of the raster
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
save(Gamma_shannon_slope_res_300, file="~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_300.Rdata")
get(load("~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_300.Rdata"))

observed <- data$Shannon_index
fit <- Gamma_shannon_slope_res_300$summary.fitted.values$mean[1:length(observed)]
pseudo_r2=function(observed,fit){
  res =  observed-fit
  RRes=sum((res)^2,na.rm = T)
  RRtot=sum((observed-mean(fit,na.rm=T))^2,na.rm = T)
  pseudo_r2_val=1-RRes/RRtot
  print(RRes)
  print(RRtot)
  return(pseudo_r2_val)  
}
pr <- round(pseudo_r2(observed, fit), digits = 3)
qq300 <- ggplot()+  
  geom_abline(slope=1, intercept=0, col="grey")+
  geom_point(data=data.frame(observed=observed, fit=fit), aes(fit, observed), pch=21)+
  ylab("")+
  xlab("")+
  scale_x_continuous(breaks=c(1,2,3), limits =c(1, 3.5))+
  scale_y_continuous(breaks=c(1,2,3), limits =c(1, 4))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
        panel.grid.major.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=margin(t = 1, r = 100, b = 1, l = 1, unit = "pt"), 
        panel.background = element_rect(fill=NA),
        text = element_text(size=20))+
  draw_text(text=paste("pseudo-R2 =", pr), x=1.7, y=3, size=20)
qq300
-sum(log(Gamma_shannon_slope_res_300$cpo$cpo))

qqplot_sd_slope_together <- plot_grid(qq100, qq200, qq300, align = "v", ncol=1, axis = "r", labels="AUTO", label_size = 25)
# qqplot_sd_slope_together
# save_plot(qqplot_sd_slope_together, filename = "~/data/output/final_plot/qqplot_sd_slope_together_aviris.png", base_height = 25, base_width = 10,bg="white")

library(grid)
library(gridExtra)
x.grob <- textGrob(paste("Fitted Shannon index from AVIRIS-NG data, modelled", "with standard variation of topographic slope as predictor", sep="\n"), gp=gpar(fontsize=25), vjust = 0.2)
y.grob <- textGrob("Observe Shannon index from AVIRIS-NG data", rot=90, gp=gpar(fontsize=25)) 
qqplot_sd_slope_together_final <- grid.arrange(arrangeGrob(qqplot_sd_slope_together,bottom = x.grob, left = y.grob))
qqplot_sd_slope_together_final
save_plot(qqplot_sd_slope_together_final, filename = "~/data/output/final_plot/qqplot_sd_slope_together_final_aviris.png", base_height = 25, base_width = 11,bg="white")

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
for(x in 1:length(fact)){
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


###########
# Modelling sd(slope)
###########

# 100m 
window_size <- c(10, 20, 30)
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
Gamma_shannon_slope_res_100_sent <- inla(formula,
                                    family = "gamma",
                                    data = data,
                                    control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_slope_res_100_sent)
save(Gamma_shannon_slope_res_100_sent, file="~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_100_sent.Rdata")
get(load("~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_100_sent.Rdata"))

observed <- data$Shannon_index
fit <- Gamma_shannon_slope_res_100_sent$summary.fitted.values$mean[1:length(observed)]
pseudo_r2=function(observed,fit){
  res =  observed-fit
  RRes=sum((res)^2,na.rm = T)
  RRtot=sum((observed-mean(fit,na.rm=T))^2,na.rm = T)
  pseudo_r2_val=1-RRes/RRtot
  print(RRes)
  print(RRtot)
  return(pseudo_r2_val)  
}
pr <- round(pseudo_r2(observed, fit), digits = 3)
qq100 <- ggplot()+  
  geom_abline(slope=1, intercept=0, col="grey")+
  geom_point(data=data.frame(observed=observed, fit=fit), aes(fit, observed), pch=21)+
  ylab("")+
  xlab("")+
  scale_x_continuous(breaks=c(1,2,3), limits =c(1, 3.5))+
  scale_y_continuous(breaks=c(1,2,3), limits =c(1, 4))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
        panel.grid.major.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=margin(t = 1, r = 4, b = 1, l = 1, unit = "pt"), 
        panel.background = element_rect(fill=NA),
        text = element_text(size=20))+
  draw_text(text=paste("pseudo-R2 =", pr), x=1.7, y=3, size=20)
qq100
-sum(log(Gamma_shannon_slope_res_100_sent$cpo$cpo))


# 200m 
window_size <- c(10, 20, 30)
x <- 2
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_", window_size[x], "_PC1278", sep=""))
plot(Shannon)
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast(paste("~/data/ArcDEM/sd_slope_sent_aviris_aoi_", window_size[x], "_res.tif", sep=""))
Slope_matrix <- as.matrix(Slope, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
length(which(s==0))
s[!is.na(s)] <- s[!is.na(s)]+ 1
sl <- inla.matrix2vector(Slope_matrix)
table(is.na(Slope_matrix))
table(is.na(Shannon_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, slope = sl, node=node)
data <- data[-c(which(is.na(s))),]
hist(data$Shannon_index)
plot(data$Shannon_index~log(data$slope))

formula= Shannon_index ~ 1 + slope +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_slope_res_200_sent <- inla(formula,
                                    family = "gamma",
                                    data = data,
                                    control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_slope_res_200_sent)
save(Gamma_shannon_slope_res_200_sent, file="~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_200_sent.Rdata")
get(load("~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_200_sent.Rdata"))

observed <- data$Shannon_index
fit <- Gamma_shannon_slope_res_200_sent$summary.fitted.values$mean[1:length(observed)]
pseudo_r2=function(observed,fit){
  res =  observed-fit
  RRes=sum((res)^2,na.rm = T)
  RRtot=sum((observed-mean(fit,na.rm=T))^2,na.rm = T)
  pseudo_r2_val=1-RRes/RRtot
  print(RRes)
  print(RRtot)
  return(pseudo_r2_val)  
}
pr <- round(pseudo_r2(observed, fit), digits = 3)
qq200 <- ggplot()+  
  geom_abline(slope=1, intercept=0, col="grey")+
  geom_point(data=data.frame(observed=observed, fit=fit), aes(fit, observed), pch=21)+
  ylab("")+
  xlab("")+
  scale_x_continuous(breaks=c(1,2,3), limits =c(1, 3.5))+
  scale_y_continuous(breaks=c(1,2,3), limits =c(1, 4))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
        panel.grid.major.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=margin(t = 1, r = 4, b = 1, l = 1, unit = "pt"), 
        panel.background = element_rect(fill=NA),
        text = element_text(size=20))+
  draw_text(text=paste("pseudo-R2 =", pr), x=1.7, y=3, size=20)
qq200
-sum(log(Gamma_shannon_slope_res_200_sent$cpo$cpo))

# 300m 
window_size <- c(10, 20, 30)
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
Gamma_shannon_slope_res_300_sent <- inla(formula,
                                    family = "gamma",
                                    data = data,
                                    control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_slope_res_300_sent)
save(Gamma_shannon_slope_res_300_sent, file="~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_300_sent.Rdata")
get(load("~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_300_sent.Rdata"))

observed <- data$Shannon_index
fit <- Gamma_shannon_slope_res_300_sent$summary.fitted.values$mean[1:length(observed)]
pseudo_r2=function(observed,fit){
  res =  observed-fit
  RRes=sum((res)^2,na.rm = T)
  RRtot=sum((observed-mean(fit,na.rm=T))^2,na.rm = T)
  pseudo_r2_val=1-RRes/RRtot
  print(RRes)
  print(RRtot)
  return(pseudo_r2_val)  
}
pr <- round(pseudo_r2(observed, fit), digits = 3)
qq300 <- ggplot()+  
  geom_abline(slope=1, intercept=0, col="grey")+
  geom_point(data=data.frame(observed=observed, fit=fit), aes(fit, observed), pch=21)+
  ylab("")+
  xlab("")+
  scale_x_continuous(breaks=c(1,2,3), limits =c(1, 3.5))+
  scale_y_continuous(breaks=c(1,2,3), limits =c(1, 4))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
        panel.grid.major.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=margin(t = 1, r = 50, b = 2, l = 1, unit = "pt"), 
        panel.background = element_rect(fill=NA),
        text = element_text(size=20))+
  draw_text(text=paste("pseudo-R2 =", pr), x=1.7, y=3, size=20)
qq300
-sum(log(Gamma_shannon_slope_res_300_sent$cpo$cpo))


qqplot_sd_slope_together <- plot_grid(qq100, qq200, qq300, align = "v", ncol=1, axis = "r", labels="AUTO", label_size = 25)
# qqplot_sd_slope_together
# save_plot(qqplot_sd_slope_together, filename = "~/data/output/final_plot/qqplot_sd_slope_together_sent.png", base_height = 25, base_width = 10,bg="white")

library(grid)
library(gridExtra)
x.grob <- textGrob(paste("Fitted Shannon index from Sentinel-2 data, modelled", "with standard variation of topographic slope as predictor", sep="\n"), gp=gpar(fontsize=25),vjust = 0.2)
y.grob <- textGrob("Observe Shannon index from AVIRIS-NG data", rot=90, gp=gpar(fontsize=25)) 
qqplot_sd_slope_together_sent <- grid.arrange(arrangeGrob(qqplot_sd_slope_together,bottom = x.grob, left = y.grob))
qqplot_sd_slope_together_sent
save_plot(qqplot_sd_slope_together_sent, filename = "~/data/output/final_plot/qqplot_sd_slope_together_final_sent.png", base_height = 25, base_width = 11,bg="white")


############
# plotting model coefficent estimates 
############

### coefficient of model derived from AVIRIS-NG data

get(load("~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_100.Rdata"))
get(load("~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_200.Rdata"))
get(load("~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_300.Rdata"))

Slope_df_aviris <- rbind(Gamma_shannon_slope_res_100$summary.fixed[2,c(1,3,5)],
                  Gamma_shannon_slope_res_200$summary.fixed[2,c(1,3,5)],
                  Gamma_shannon_slope_res_300$summary.fixed[2,c(1,3,5)])
est <- c("100 m", "200 m", "300 m")
Slope_df_aviris <- cbind(as.factor(est), Slope_df_aviris)
colnames(Slope_df_aviris)[c(1,3:4)] <- c("resolution","lower", "upper")
# Slope_df_aviris[,2:4] <- exp(Slope_df_aviris[,2:4])
Slope_df_aviris$resolution <- factor(Slope_df_aviris$resolution, levels = c("100 m", "200 m", "300 m"))

sd_slope_aviris_plot <- ggplot(data = Slope_df_aviris, 
                              aes(x = mean, y = resolution, xmin = lower, xmax = upper)) +
  geom_pointrange(fatten=2.5) +
  labs(
    # title = expression(paste("Aviris data: Model estimate of ", sigma, "(slope) on the estimate Shannon index")),
       x = expression(atop(paste("Coefficient estimates of ", sigma, " (slope) to predict"), 
                           " AVIRIS-NG-derived spectral diverstiy", sep="\n")),
       y = ""
       # caption = "Models fit with INLA. Error bars show the 95% confidence interval."
       )+
  geom_vline(xintercept = as.numeric(0), col="red", alpha=0.5, lty="dashed")+
  theme_bw() + 
  coord_cartesian(xlim = c(-0.025, 0.065))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.y = element_text(size=12), axis.text.x = element_text(size=12)
        ,panel.grid.major.x = element_blank(),  panel.grid.major.y = element_blank(), 
        axis.title.x = element_text(vjust=-1.8, hjust=0.3, size=14), 
        plot.margin=margin(t = 5.5, r = 8, b = 6, l = 1, unit = "pt"), 
  )
sd_slope_aviris_plot


### coefficient of model derived from Sentinel-2 data

get(load("~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_100_sent.Rdata"))
get(load("~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_200_sent.Rdata"))
get(load("~/scratch/INLA_modelling/Aviris_modelling/Gamma_shannon_slope_res_300_sent.Rdata"))

Slope_df_sent <- rbind(Gamma_shannon_slope_res_100_sent$summary.fixed[2,c(1,3,5)],
                         Gamma_shannon_slope_res_200_sent$summary.fixed[2,c(1,3,5)],
                         Gamma_shannon_slope_res_300_sent$summary.fixed[2,c(1,3,5)])
est <- c("100 m", "200 m", "300 m")
Slope_df_sent <- cbind(as.factor(est), Slope_df_sent)
colnames(Slope_df_sent)[c(1,3:4)] <- c("resolution","lower", "upper")
# Slope_df_sent[,2:4] <- exp(Slope_df_sent[,2:4])
Slope_df_sent$resolution <- factor(Slope_df_sent$resolution, levels = c("100 m", "200 m", "300 m"))

Slope_df_sent_plot <- ggplot(data = Slope_df_sent, 
                              aes(x = mean, y = resolution, xmin = lower, xmax = upper)) +
  geom_pointrange(fatten=2.5) +
  labs(
    # title = expression(paste("Sentinel-2 data: Model estimate of ", sigma, "(slope) on the estimate Shannon index")),
       x = expression(atop(paste("Coefficient estimates of ", sigma, " (slope) to predict"), 
                           " Sentinel-2-derived spectral diverstiy", sep="\n")),
       y = ""
       # caption = "Models fit with INLA. Error bars show the 95% confidence interval."
       )+
  geom_vline(xintercept = as.numeric(0), col="red", alpha=0.5, lty="dashed")+
  theme_bw() + 
  coord_cartesian(xlim = c(-0.025, 0.065))+
  # annotate("text", x=1, y=1, label= "pseudo-R2", cex=0.5, col="red") + 
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.y = element_text(size=12), axis.text.x = element_text(size=12)
        ,panel.grid.major.x = element_blank(),  panel.grid.major.y = element_blank(), 
        axis.title.x = element_text(vjust=-1.8, hjust=0.3, size=14), 
        plot.margin=margin(t = 5.5, r = 8, b = 6, l = 1, unit = "pt")
  )
Slope_df_sent_plot

grob <- grobTree(textGrob("pseudo-R2", x=4,  y=0.25, hjust=0.5,
                          gp=gpar(col="red", fontsize=13, fontface="italic")))
Slope_df_sent_plot + annotation_custom(grob)


### Putting plots together

sd_slope_aviris <- plot_grid(sd_slope_aviris_plot, Slope_df_sent_plot)
x.grob <- textGrob("")
y.grob <- textGrob("Resolution", rot=90, gp=gpar(fontsize=14)) 
final_plot <- grid.arrange(arrangeGrob(sd_slope_aviris,bottom = x.grob, left = y.grob))
save_plot(final_plot, filename="~/data/output/final_plot/coeff_plot_sd_slope_aviris.png", ncol=2, base_height = 4, base_width = 5)


# ######
# # semi-variogram analysis
# ######
# library(gstat)
# library(sp)
# library(raster)
# aviris_map <- rast("~/data/biodivmapR_Aviris/RESULTS/ang20190802t220708_rfl_rect/SPCA/ALPHA/Shannon_20_PC1289")
# aviris_map_raster <- raster(aviris_map)
# aviris_map_sp <- as(aviris_map_raster, "SpatialPointsDataFrame")
# vario_close <- variogram(Shannon_20_PC1289~1, data=aviris_map_sp, cutoff=5000, width=100)
# plot(vario_close)
# try <- vgm(0.1, "Exp", 1500, nugget=0.05)
# plot(try, cutoff=5000)
# vario_fit_close <- fit.variogram(vario_close, try, fit.kappa = TRUE)
# plot(vario_fit_close, cutoff=5000)
# plot(vario_close, vario_fit_close)
# vario_fit_close
# ggplot()+
#   geom_point(data=vario_close, aes(dist, gamma))+
#   theme_cowplot()
# save_plot(plot, filename = "~/data/output/final_plot/variogram_close_shannon_res100_width100_cutoff5000.png", base_height = 4)










