setwd("~/data")
library(sf)
library(tidyverse)
library(dplyr)
library(tidyterra)
library(remotes)
library(dissUtils)
library(ggplot2)
library(cowplot)

library(sp)
library(rgdal)  
require(spdep)
library(INLA)
library(terra)
library(INLAutils)
library(scales)

#################
# Scaling analysis 
#################

# Get Shannon index data at different scale

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
# Apply normalization with continuum removal?
Continuum_Removal <- F
# Type of dimensionality reduction
TypePCA <- 'SPCA'
# PCA FILTERING:        Set to TRUE if you want second filtering based on PCA outliers to be processed.
# Slower process
# Automatically set to FALSE if TypePCA     = 'MNF'
FilterPCA <- F
nbCPU <- 2
MaxRAM <- 12
# number of clusters (spectral species)
nbclusters <- 20

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
# assessment

viridis_colors <- viridis::plasma(20)
for(x in 1:length(window_size)){
  path <- paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_", window_size[x],"_PC1278", sep="")
  m <- rast(path)
  plot(m, col=viridis_colors)
    
}


# prepare elevation data at different scale

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
  tile_DEM_crop_aggregated <- aggregate(tile_DEM_crop_resample_noNA, fact=window_size[x], fun="sd")
  tile_DEM_masked <- mask(tile_DEM_crop_aggregated, Shannon)
  writeRaster(tile_DEM_masked, filename = paste("~/data/ArcDEM/ArcDEM_masked_", window_size[x], "_res.tif", sep=""))
}


#################
# Modelling sd(elevation)
#################

# resolution 100m x 100m
# Sentinel_lattice <-readOGR("~/data/output/INLA_modelling/model_object_res_10_with_NA.shp")
# Sentinel_data <- Sentinel_lattice@data
# Sentinel_data$Shannon_10[!is.na(Sentinel_data$Shannon_10)] <- Sentinel_data$Shannon_10[!is.na(Sentinel_data$Shannon_10)]+ 1
# which(Sentinel_data$Shannon_10==0)
# hist(Sentinel_data$Shannon_10)
# hist(Sentinel_data$sd_topo)
# hist(log(Sentinel_data$sd_topo))
# plot(Sentinel_data$Shannon_10~Sentinel_data$sd_topo)
# plot(Sentinel_data$Shannon_10~log(Sentinel_data$sd_topo))

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
data <- data.frame(Shannon_index = s, topo = e, node=node)
data <- data[-c(which(is.na(s))),]

# log.range = list(initial = log(5), fixed=TRUE) 
# hyperpar_matern = list(initial = -3, param=c(23.36,0.001))
# formula_matern = y ~ wildfire + logging + lichwood + openlich + deciduous + 
#   water + wetland + meanelev +
#   f(node_matern, model = "matern2d", nrow = nrow.larger, 
#     ncol = ncol.larger, hyper = list(range = log.range, prec 
#                                      = hyperpar_matern))
# hyper = list(range = list(param =c(1, 1),prior = "loggamma",initial=1),
#              prec = list(param=c(1, 1)))

formula= Shannon_index ~ 1+ log(topo)+
  f(node, model="matern2d", nrow=nrow, ncol=ncol)

Gaussian_shannon_log_topo <- inla(formula,family = "gaussian",
                                  control.family=list(link='identity'),
                                  data = data,
                                  control.compute = list(cpo = T, dic = T, waic = T, return.marginals.predictor=TRUE, config=T), verbose=TRUE)
summary(Gaussian_shannon_log_topo)
save(Gaussian_shannon_log_topo, file="~/scratch/INLA_modelling/Gaussian_shannon_log_topo.Rdata")
# Gaussian_loglink <- get(load("~/scratch/INLA_modelling/Gaussian_model_loglink_shannon_log_topo.Rdata"))
observed=data$Shannon_index
plot_inla_residuals(Gaussian_shannon_log_topo, observed = observed)
# 
# 
# plot_inla_residuals(Gaussian_loglink_withNA, observed=data$Shannon_index)
# 
# ggplot_inla_residuals(Model_Lattice_no_ndvi, observed=observed)
# ggplot_inla_residuals2(Model_Lattice_no_ndvi, observed, se = FALSE)
# 
# 
# formula= Shannon_index ~ 1+ log(topo)+
#   f(node, model="matern2d", nrow=nrow, ncol=ncol)
# 
# Gaussian_shannon_log_topo <- inla(formula,family = "gaussian",
#                                       control.family=list(link='identity'),
#                                       data = data,
#                                       control.compute = list(cpo = T, dic = T, waic = T, return.marginals.predictor=TRUE), verbose=TRUE)
# 
# 
# summary(Gaussian_shannon_log_topo)
# plot(data$Shannon_index~log(data$topo))

#####
# Gamma
#####
# seems better than gaussian

formula= Shannon_index ~ 1+ topo+
  f(node, model="matern2d", nrow=nrow, ncol=ncol)

Gamma_shannon_topo <- inla(formula,     
                               family = "gamma",
                               data = data,
                               control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)

summary(Gamma_shannon_topo)

Gamma_shannon_topo_res100 <- Gamma_shannon_topo
save(Gamma_shannon_topo_res100, file="~/data/output/INLA_modelling/Gamma_shannon_topo_res100.Rdata")
get(load("~/scratch/INLA_modelling/Gamma_shannon_topo_res100.Rdata"))

observed <- data$Shannon_index
plot_inla_residuals_ML(Gamma_shannon_topo_res100, observed=observed)

#################
# resolution 200mx200m
#################
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 20, "PC1278", sep="_"))
Elev <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 20, "_res.tif", sep=""))
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
data <- data.frame(Shannon_index = s, topo = e, node=node)
data <- data[-c(which(is.na(s))),]


formula= Shannon_index ~ 1+ topo+
  f(node, model="matern2d", nrow=nrow, ncol=ncol)

Gamma_shannon_topo_res200 <- inla(formula,     
                           family = "gamma",
                           data = data,
                           control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)

summary(Gamma_shannon_topo_res200)

save(Gamma_shannon_topo_res200, file="~/data/output/INLA_modelling/Gamma_shannon_topo_res200.Rdata")
get(load("~/scratch/INLA_modelling/Gamma_shannon_topo_res200.Rdata"))

observed <- log(data$Shannon_index) 
plot_inla_residuals(Gamma_shannon_topo_res200, observed=observed)

###############
# resolution 300m x 300m
###############
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 30, "PC1278", sep="_"))
Elev <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 30, "_res.tif", sep=""))
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
data <- data.frame(Shannon_index = s, topo = e, node=node)
data <- data[-c(which(is.na(s))),]


formula= Shannon_index ~ 1+ topo+
  f(node, model="matern2d", nrow=nrow, ncol=ncol)

Gamma_shannon_topo_res300 <- inla(formula,     
                                  family = "gamma",
                                  data = data,
                                  control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)

summary(Gamma_shannon_topo_res300)

save(Gamma_shannon_topo_res300, file="~/scratch/INLA_modelling/Gamma_shannon_topo_res300.Rdata")
get(load("~/scratch/INLA_modelling/Gamma_shannon_topo_res300.Rdata"))

observed <- log(data$Shannon_index) 
plot_inla_residuals(Gamma_shannon_topo_res300, observed=observed)

################
# resolution 500m x 500m
################
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 50, "PC1278", sep="_"))
Elev <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 50, "_res.tif", sep=""))
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
data <- data.frame(Shannon_index = s, topo = e, node=node)
data <- data[-c(which(is.na(s))),]


formula= Shannon_index ~ 1+ topo+
  f(node, model="matern2d", nrow=nrow, ncol=ncol)

Gamma_shannon_topo_res500 <- inla(formula,     
                                  family = "gamma",
                                  data = data,
                                  control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)

summary(Gamma_shannon_topo_res500)

save(Gamma_shannon_topo_res500, file="~/data/output/INLA_modelling/Gamma_shannon_topo_res500.Rdata")
get(load("~/scratch/INLA_modelling/Gamma_shannon_topo_res500.Rdata"))


##############
# resolution 1000m x 1000m
##############
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 100, "PC1278", sep="_"))
Elev <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 100, "_res.tif", sep=""))
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
data <- data.frame(Shannon_index = s, topo = e, node=node)
data <- data[-c(which(is.na(s))),]


formula= Shannon_index ~ 1+ topo+
  f(node, model="matern2d", nrow=nrow, ncol=ncol)

Gamma_shannon_topo_res1000 <- inla(formula,     
                                  family = "gamma",
                                  data = data,
                                  control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)

summary(Gamma_shannon_topo_res1000)
Gamma_shannon_topo_res1000_sum <- summary(Gamma_shannon_topo_res1000)
Gamma_shannon_topo_res_10000_df  <- data.frame(mean = Gamma_shannon_topo_res1000_sum$fixed[2,"mean"],
                       lower = Gamma_shannon_topo_res1000_sum$fixed[2,"0.025quant"],
                       upper = Gamma_shannon_topo_res1000_sum$fixed[2,"0.975quant"],
                       stringsAsFactors = FALSE)

est <- c("sd(topography)")
Gamma_shannon_topo_res_10000_df <- bind_cols(as.factor(est), Gamma_shannon_topo_res_10000_df)

sd_topo_res_1000_plot <- Gamma_shannon_topo_res_10000_df %>% 
  mutate(est= factor(est, levels= c("sd(topography)"))) %>% 
  ggplot() + 
  geom_pointrange(aes(x = factor(est), y = mean, ymin = lower, ymax = upper), position = position_dodge(0.5)) +
  xlab("coefficient") +
  ylab("Estimate")+
  coord_flip()+
  theme_bw()

sd_topo_res_1000_plot 

save(Gamma_shannon_topo_res1000, file="~/data/output/INLA_modelling/Gamma_shannon_topo_res1000.Rdata")
get(load("~/scratch/INLA_modelling/Gamma_shannon_topo_res1000.Rdata"))

#################
# plotting
#################
Gamma_shannon_topo_res1000_sum <- summary(Gamma_shannon_topo_res1000)
Gamma_shannon_topo_res1000_df  <- data.frame(mean = Gamma_shannon_topo_res1000_sum$fixed[2,"mean"],
                                               lower = Gamma_shannon_topo_res1000_sum$fixed[2,"0.025quant"],
                                               upper = Gamma_shannon_topo_res1000_sum$fixed[2,"0.975quant"],
                                               stringsAsFactors = FALSE)

Gamma_shannon_topo_res100_sum <- summary(Gamma_shannon_topo_res100)
Gamma_shannon_topo_res100_df  <- data.frame(mean = Gamma_shannon_topo_res100_sum$fixed[2,"mean"],
                                               lower = Gamma_shannon_topo_res100_sum$fixed[2,"0.025quant"],
                                               upper = Gamma_shannon_topo_res100_sum$fixed[2,"0.975quant"],
                                               stringsAsFactors = FALSE)

Gamma_shannon_topo_res200_sum <- summary(Gamma_shannon_topo_res200)
Gamma_shannon_topo_res200_df  <- data.frame(mean = Gamma_shannon_topo_res200_sum$fixed[2,"mean"],
                                            lower = Gamma_shannon_topo_res200_sum$fixed[2,"0.025quant"],
                                            upper = Gamma_shannon_topo_res200_sum$fixed[2,"0.975quant"],
                                            stringsAsFactors = FALSE)

Gamma_shannon_topo_res300_sum <- summary(Gamma_shannon_topo_res300)
Gamma_shannon_topo_res300_df  <- data.frame(mean = Gamma_shannon_topo_res300_sum$fixed[2,"mean"],
                                            lower = Gamma_shannon_topo_res300_sum$fixed[2,"0.025quant"],
                                            upper = Gamma_shannon_topo_res300_sum$fixed[2,"0.975quant"],
                                            stringsAsFactors = FALSE)

Gamma_shannon_topo_res500_sum <- summary(Gamma_shannon_topo_res500)
Gamma_shannon_topo_res500_df  <- data.frame(mean = Gamma_shannon_topo_res500_sum$fixed[2,"mean"],
                                            lower = Gamma_shannon_topo_res500_sum$fixed[2,"0.025quant"],
                                            upper = Gamma_shannon_topo_res500_sum$fixed[2,"0.975quant"],
                                            stringsAsFactors = FALSE)

est <- c("100m resolution", "200m resolution", "300m resolution", "500m resolution", "1000m resolution")
Gamma_shannon_topo_df <- rbind(Gamma_shannon_topo_res100_df, Gamma_shannon_topo_res200_df, Gamma_shannon_topo_res300_df, Gamma_shannon_topo_res500_df, Gamma_shannon_topo_res1000_df)
Gamma_shannon_topo_df <- bind_cols(as.factor(est), Gamma_shannon_topo_df)

Gamma_shannon_topo_df_trans <- Gamma_shannon_topo_df
Gamma_shannon_topo_df_trans[,2:4] <- exp(Gamma_shannon_topo_df_trans[2:4])

sd_topo_plot <- Gamma_shannon_topo_df %>% 
  mutate(est= factor(est, levels= c("100m resolution", "200m resolution", "300m resolution", "500m resolution", "1000m resolution"))) %>% 
  ggplot() + 
  geom_pointrange(aes(x = factor(est), y = mean, ymin = lower, ymax = upper), position = position_dodge(0.5)) +
  xlab("coefficient sd(elevation)") +
  ylab("Estimate")+
  coord_flip()+
  geom_vline(xintercept = 0, linetype='longdash', col="red")+
  theme_bw()

sd_topo_plot

sd_topo_plot_trans <- Gamma_shannon_topo_df_trans %>% 
  mutate(est= factor(est, levels= c("100m resolution", "200m resolution", "300m resolution", "500m resolution", "1000m resolution"))) %>% 
  ggplot() + 
  geom_pointrange(aes(x = factor(est), y = mean, ymin = lower, ymax = upper), position = position_dodge(0.5)) +
  labs(title = expression(paste("Model estimates of ", sigma, "(elevation) on the estimate Shannon index")),
       x = "Resolution",
       y = "Coefficient Estimate",
       caption = "Models fit with INLA. Error bars show the 95% confidence interval.")+
  coord_flip()+
  theme_bw()+
  geom_hline(yintercept = as.numeric(1), linetype='longdash', col="red", alpha=0.5)
sd_topo_plot_trans
geom_vline(data=ev,aes(xintercept=as.numeric(dt)))

par(mfrow=c(1,2))
sd_topo_plot
sd_topo_plot_trans
plot_grid(sd_topo_plot, sd_topo_plot_trans)

library(sjPlot)
plot_model(Gamma_shannon_topo_res100)


#################
# sd(ele) with 10m aggregated data
#################
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 10, "PC1278", sep="_"))
Elev <- rast(paste("~/data/ArcDEM/sd_ele_agg10_", 10, "_res.tif", sep=""))
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
data <- data.frame(Shannon_index = s, topo = e, node=node)
data <- data[-c(which(is.na(s))),]

formula= Shannon_index ~ 1+ topo+
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_topo_agg10_res100 <- inla(formula,     
                                  family = "gamma",
                                  data = data,
                                  control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)

summary(Gamma_shannon_topo_agg10_res100)
save(Gamma_shannon_topo_agg10_res100, file="~/data/output/INLA_modelling/Gamma_shannon_topo_agg10_res100.Rdata")
get(load("~/data/output/INLA_modelling/Gamma_shannon_topo_agg10_res100.Rdata"))

# 200m
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 20, "PC1278", sep="_"))
Elev <- rast(paste("~/data/ArcDEM/sd_ele_agg10_", 20, "_res.tif", sep=""))
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
data <- data.frame(Shannon_index = s, topo = e, node=node)
data <- data[-c(which(is.na(s))),]

formula= Shannon_index ~ 1+ topo+
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_topo_agg10_res200 <- inla(formula,     
                                        family = "gamma",
                                        data = data,
                                        control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)

summary(Gamma_shannon_topo_agg10_res200)
save(Gamma_shannon_topo_agg10_res200, file="~/data/output/INLA_modelling/Gamma_shannon_topo_agg10_res200.Rdata")
get(load("~/data/output/INLA_modelling/Gamma_shannon_topo_agg10_res200.Rdata"))


# 300m
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 30, "PC1278", sep="_"))
Elev <- rast(paste("~/data/ArcDEM/sd_ele_agg10_", 30, "_res.tif", sep=""))
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
data <- data.frame(Shannon_index = s, topo = e, node=node)
data <- data[-c(which(is.na(s))),]

formula= Shannon_index ~ 1+ topo+
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_topo_agg10_res300 <- inla(formula,     
                                        family = "gamma",
                                        data = data,
                                        control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)

summary(Gamma_shannon_topo_agg10_res300)
save(Gamma_shannon_topo_agg10_res300, file="~/data/output/INLA_modelling/Gamma_shannon_topo_agg10_res300.Rdata")
get(load("~/data/output/INLA_modelling/Gamma_shannon_topo_agg10_res200.Rdata"))

# 1000m
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 100, "PC1278", sep="_"))
Elev <- rast(paste("~/data/ArcDEM/sd_ele_agg10_", 100, "_res.tif", sep=""))
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
data <- data.frame(Shannon_index = s, topo = e, node=node)
data <- data[-c(which(is.na(s))),]

formula= Shannon_index ~ 1+ topo+
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_topo_agg10_res1000 <- inla(formula,     
                                        family = "gamma",
                                        data = data,
                                        control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)

summary(Gamma_shannon_topo_agg10_res1000)
save(Gamma_shannon_topo_agg10_res1000, file="~/data/output/INLA_modelling/Gamma_shannon_topo_agg10_res1000.Rdata")
get(load("~/data/output/INLA_modelling/Gamma_shannon_topo_agg10_res500.Rdata"))

# 500m
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 50, "PC1278", sep="_"))
Elev <- rast(paste("~/data/ArcDEM/sd_ele_agg10_", 50, "_res.tif", sep=""))
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
data <- data.frame(Shannon_index = s, topo = e, node=node)
data <- data[-c(which(is.na(s))),]

formula= Shannon_index ~ 1+ topo+
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_topo_agg10_res500 <- inla(formula,     
                                        family = "gamma",
                                        data = data,
                                        control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)

summary(Gamma_shannon_topo_agg10_res500)
save(Gamma_shannon_topo_agg10_res500, file="~/data/output/INLA_modelling/Gamma_shannon_topo_agg10_res500.Rdata")
get(load("~/data/output/INLA_modelling/Gamma_shannon_topo_agg10_res200.Rdata"))


##################
# simple glm regression for sd(elevation)
##################
# Generate some negative log data with some gaussian noise
# 3 samples per data point
set.seed(10)
my_data <- data.frame(
  x = 1:100,
  y = -log(1:100) + rnorm(n = 300)
)

# Check data summary
summary(my_data)
hist(my_data$y)

# Y has negative values so push up by 10 to allow for gaussian log link glm
my_data$y <- my_data$y + 10

# Plot of data
(my_plot <- ggplot(
  my_data,
  aes(x = x, y = y)) +
    geom_point() +
    scale_y_continuous(limits = c(0, 15)) +
    theme_cowplot())


# Fit log link gaussian glm
my_glm <- glm(y ~ x, family = gaussian(link = "log"), data = my_data)
summary(my_glm)

# Get predictions from inbuild function for x = 1:100
my_preds <- predict(
  my_glm,
  newdata = data.frame(x = 1:100),
  type = "response"
)

# Calculate predictions ourselves:
# log(y) = a * x + b
# => y = exp(a * x + b)
coef(my_glm)
a <- coef(my_glm)[2]
b <- coef(my_glm)[1]
y_preds <- exp(a * 1:100 + b)
names(y_preds) <- 1:100

# compare
identical(my_preds, y_preds)

# generate predctions data frame
my_preds_df <- data.frame(
  x = as.numeric(names(my_preds)),
  y = my_preds
)

# Plot preds and true relationship
my_plot +
  geom_line(data = my_preds_df, 
            colour = "blue")


# for 100m resolution
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
data <- data.frame(Shannon_index = s, topo = e, node=node)
data <- data[-c(which(is.na(s))),]

res_100_glm <- glm(Shannon_index ~ topo, family = Gamma(link = "log"), data = data)
summary(res_100_glm)
plot(data$Shannon_index~data$topo)
(res_100_plot <- ggplot(
  data,
  aes(x = topo, y = Shannon_index)) +
    geom_point(shape=21) +
    scale_y_continuous(limits = c(1, 4)) +
    labs(x = "sd(elevation)", y="Shannon index", title = "100m resolution")+
    theme_cowplot())

res_100_pred <- predict(
  res_100_glm,
  newdata = data.frame(topo = seq(0, 20, 0.01)),
  type = "response"
)


res_100_pred_df <- data.frame(
  topo = seq(0,20,0.01),
  Shannon_index = res_100_pred
)

# Plot preds and true relationship
res_100_plot_pred <- res_100_plot +
  geom_line(data = res_100_pred_df, 
            colour = "blue")
res_100_plot_pred

# # Gamma regression with identity link
# data$log_topo <- log(data$topo)
# res_100_glm <- glm(Shannon_index ~ log_topo, family = Gamma(link = "identity"), data = data)
# summary(res_100_glm)
# (res_100_plot <- ggplot(
#   data,
#   aes(x = log_topo, y = Shannon_index)) +
#     geom_point(shape=21) +
#     scale_y_continuous(limits = c(1, 4)) +
#     theme_cowplot())
# 
# res_100_pred <- predict(
#   res_100_glm,
#   newdata = data.frame(log_topo = log(seq(0.01, 20, 0.01))),
#   type = "response"
# )
# res_100_pred_df <- data.frame(
#   log_topo = log(seq(0.01,20,0.01)),
#   Shannon_index = res_100_pred
# )
# res_100_plot +
#   geom_line(data = res_100_pred_df, 
#             colour = "blue")


# for 200m resolution
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 20, "PC1278", sep="_"))
Elev <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 20, "_res.tif", sep=""))
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
data <- data.frame(Shannon_index = s, topo = e, node=node)
data <- data[-c(which(is.na(s))),]

res_200_glm <- glm(Shannon_index ~ topo, family = Gamma(link = "log"), data = data)
summary(res_200_glm)
(res_200_plot <- ggplot(
  data,
  aes(x = topo, y = Shannon_index)) +
    geom_point(shape=21) +
    scale_y_continuous(limits = c(1, 4)) +
    labs(x = "sd(elevation)", y="Shannon index", title = "200m resolution")+
    theme_cowplot())
range(data$topo)*1/1000
res_200_pred <- predict(
  res_200_glm,
  newdata = data.frame(topo = seq(0, 30, 0.03)),
  type = "response"
)
res_200_pred_df <- data.frame(
  topo = seq(0,30,0.03),
  Shannon_index = res_200_pred
)

# Plot preds and true relationship
res_200_plot_pred <- res_200_plot +
  geom_line(data = res_200_pred_df, 
            colour = "blue")
res_200_plot_pred

#for 300m 
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 30, "PC1278", sep="_"))
Elev <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 30, "_res.tif", sep=""))
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
data <- data.frame(Shannon_index = s, topo = e, node=node)
data <- data[-c(which(is.na(s))),]

res_300_glm <- glm(Shannon_index ~ topo, family = Gamma(link = "log"), data = data)
summary(res_300_glm)
(res_300_plot <- ggplot(
  data,
  aes(x = topo, y = Shannon_index)) +
    geom_point(shape=21) +
    scale_y_continuous(limits = c(1, 4)) +
    labs(x = "sd(elevation)", y="Shannon index", title = "300m resolution")+
    theme_cowplot())
range(data$topo)
res_300_pred <- predict(
  res_300_glm,
  newdata = data.frame(topo = seq(0, 40, 0.03)),
  type = "response"
)
res_300_pred_df <- data.frame(
  topo = seq(0,40,0.03),
  Shannon_index = res_300_pred
)

# Plot preds and true relationship
res_300_plot_pred <-res_300_plot +
  geom_line(data = res_300_pred_df, 
            colour = "blue")
res_300_plot_pred

# for 500m
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 50, "PC1278", sep="_"))
Elev <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 50, "_res.tif", sep=""))
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
data <- data.frame(Shannon_index = s, topo = e, node=node)
data <- data[-c(which(is.na(s))),]

res_500_glm <- glm(Shannon_index ~ topo, family = Gamma(link = "log"), data = data)
summary(res_500_glm)
(res_500_plot <- ggplot(
  data,
  aes(x = topo, y = Shannon_index)) +
    geom_point(shape=21) +
    scale_y_continuous(limits = c(1, 4)) +
    labs(x = "sd(elevation)", y="Shannon index", title = "500m resolution")+
    theme_cowplot())
range(data$topo)
res_500_pred <- predict(
  res_500_glm,
  newdata = data.frame(topo = seq(0, 45, 0.04)),
  type = "response"
)
res_500_pred_df <- data.frame(
  topo = seq(0,45,0.04),
  Shannon_index = res_500_pred
)

# Plot preds and true relationship
res_500_plot_pred <- res_500_plot +
  geom_line(data = res_500_pred_df, 
            colour = "blue")
res_500_plot_pred


# for 1000m 
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 100, "PC1278", sep="_"))
Elev <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 100, "_res.tif", sep=""))
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
data <- data.frame(Shannon_index = s, topo = e, node=node)
data <- data[-c(which(is.na(s))),]

res_1000_glm <- glm(Shannon_index ~ topo, family = Gamma(link = "log"), data = data)
summary(res_1000_glm)
(res_1000_plot <- ggplot(
  data,
  aes(x = topo, y = Shannon_index)) +
    geom_point(shape=21) +
    scale_y_continuous(limits = c(1, 4)) +
    labs(x = "sd(elevation)", y="Shannon index", title = "1000m resolution")+
    theme_cowplot())
range(data$topo)
res_1000_pred <- predict(
  res_1000_glm,
  newdata = data.frame(topo = seq(0, 60, 0.05)),
  type = "response"
)
res_1000_pred_df <- data.frame(
  topo = seq(0,60,0.05),
  Shannon_index = res_1000_pred
)

# Plot preds and true relationship
res_1000_plot_pred <- res_1000_plot +
  geom_line(data = res_1000_pred_df, 
            colour = "blue")
res_1000_plot_pred

plot_grid(res_100_plot_pred, res_200_plot_pred, res_300_plot_pred, res_500_plot_pred, res_1000_plot_pred)

####################
# modelling with slope and aspect
####################
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 10, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Elev <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 10, "_res.tif", sep=""))
Elev_matrix <- as.matrix(Elev, wide=T)
Slope <- rast("~/data/ArcDEM/mean_slope_masked_10_res.tif")
Slope_matrix <- as.matrix(Slope, wide=T)
max_ele <- rast("~/data/ArcDEM/max_ele_masked_10_res.tif")
max_ele_matrix <- as.matrix(max_ele, wide=T)
Aspect <- rast("~/data/ArcDEM/mean_aspect_masked_10_res.tif")
Aspect_matrix <- as.matrix(Aspect, wide=T)
m <- c(0, 45, 1,
      45, 135, 2,
      135, 225, 3,
      225, 315, 4,
      315, 360, 1)
rclmat <- matrix(m, ncol=3, byrow = T)
Aspect_class <- classify(Aspect, rclmat)
plot(Aspect_class)
Aspect_class_matrix <- as.matrix(Aspect_class, wide=T)
as_c <- inla.matrix2vector(Aspect_class_matrix)
hist(as_c)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
sl <- inla.matrix2vector(Slope_matrix)
as <- inla.matrix2vector(Aspect_matrix)
e <- inla.matrix2vector(Elev_matrix)
max_e <- inla.matrix2vector(max_ele_matrix)
table(is.na(Slope_matrix))
table(is.na(Shannon_matrix))
table(is.na(Aspect_matrix))
table(is.na(Elev_matrix))
table(is.na(max_ele_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, topo=e, max_ele= max_e, slope = sl, aspect=as, aspect_class = as_c, node=node)
data <- data[-c(which(is.na(s))),]
data$aspect_class <- as.factor(data$aspect_class)
hist(data$slope)
hist(data$aspect)
hist(data$max_ele)

plot(data$Shannon_index~data$slope)
plot(data$Shannon_index~log(data$slope))
plot(data$Shannon_index~data$aspect)
plot(data$Shannon_index~data$slope, col=data$aspect_class)
plot(data$Shannon_index~data$max_ele)
plot(data$topo~data$slope)
abline(0,1,col="red")

formula= Shannon_index ~ 1+ slope+
  f(node, model="matern2d", nrow=nrow, ncol=ncol)

Gamma_shannon_slope_res_100 <- inla(formula,     
                           family = "gamma",
                           data = data,
                           control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)

summary(Gamma_shannon_slope_res_100)
summary(Gamma_shannon_topo_res100)


formula= Shannon_index ~ 1 + aspect_class +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_aspect_class_res_100 <- inla(formula,     
                                    family = "gamma",
                                    data = data,
                                    control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_aspect_class_res_100)


formula= Shannon_index ~ 1 + slope + aspect_class + slope:aspect_class+
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_aspect_class_res_100 <- inla(formula,     
                                           family = "gamma",
                                           data = data,
                                           control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)


formula= Shannon_index ~ 1 + topo + slope + aspect+
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_ele_slope_aspect_res_100 <- inla(formula,     
                                           family = "gamma",
                                           data = data,
                                           control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_ele_slope_aspect_res_100)


formula= Shannon_index ~ 1 + slope+
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_slope_res_100 <- inla(formula,
                                    family = "gamma",
                                    data = data,
                                    control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_slope_res_100)



formula= Shannon_index ~ 1 + slope + topo +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_ele_slope_res_100 <- inla(formula,
                                    family = "gamma",
                                    data = data,
                                    control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_ele_slope_res_100)


save(Gamma_shannon_ele_slope_res_100, file="~/data/output/INLA_modelling/Gamma_shannon_topo_meanslope_res100.Rdata")
get(load("~/data/output/INLA_modelling/Gamma_shannon_topo_meanslope_res100.Rdata"))

#############
# model with mean slope and topography
#############
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 20, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Elev <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 20, "_res.tif", sep=""))
Elev_matrix <- as.matrix(Elev, wide=T)
Mean_ele <- Elev <- rast(paste("~/data/ArcDEM/mean_ele_masked_", 20, "_res.tif", sep=""))
Mean_ele_matrix <- as.matrix(Mean_ele, wide=T)
Slope <- rast("~/data/ArcDEM/mean_slope_masked_20_res.tif")
Slope_matrix <- as.matrix(Slope, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
sl <- inla.matrix2vector(Slope_matrix)
e <- inla.matrix2vector(Elev_matrix)
me <- inla.matrix2vector(Mean_ele_matrix)
table(is.na(Slope_matrix))
table(is.na(Shannon_matrix))
table(is.na(Elev_matrix))
table(is.na(Mean_ele_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, topo=e, mean_ele = me, slope = sl, node=node)
data <- data[-c(which(is.na(s))),]
data$topo_by_mean <- data$topo/data$mean_ele
data$slope_by_mean <- data$slope/data$mean_ele
hist(data$slope)
hist(data$topo_by_mean)
hist(data$topo)
formula= Shannon_index ~ 1 + slope + topo +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_ele_slope_res_200 <- inla(formula,
                                        family = "gamma",
                                        data = data,
                                        control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_ele_slope_res_200)

save(Gamma_shannon_ele_slope_res_200, file="~/data/output/INLA_modelling/Gamma_shannon_ele_slope_res_200.Rdata")
get(load("~/data/output/INLA_modelling/Gamma_shannon_ele_slope_res_200.Rdata"))

glm_slope_topo <- glm(Shannon_index ~ 1 + topo+slope, family = Gamma(link = "log"), data = data)
summary(glm_slope_topo)

## 300m
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 30, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Elev <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 30, "_res.tif", sep=""))
Elev_matrix <- as.matrix(Elev, wide=T)
Slope <- rast("~/data/ArcDEM/mean_slope_masked_30_res.tif")
Slope_matrix <- as.matrix(Slope, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
sl <- inla.matrix2vector(Slope_matrix)
e <- inla.matrix2vector(Elev_matrix)
table(is.na(Slope_matrix))
table(is.na(Shannon_matrix))
table(is.na(Elev_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, topo=e, slope = sl, node=node)
data <- data[-c(which(is.na(s))),]
hist(data$slope)
formula= Shannon_index ~ 1 + slope + topo +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_ele_slope_res_300 <- inla(formula,
                                        family = "gamma",
                                        data = data,
                                        control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_ele_slope_res_300)

save(Gamma_shannon_ele_slope_res_300, file="~/data/output/INLA_modelling/Gamma_shannon_ele_slope_res_300.Rdata")
get(load("~/data/output/INLA_modelling/Gamma_shannon_ele_slope_res_300.Rdata"))

## 500m
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 50, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Elev <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 50, "_res.tif", sep=""))
Elev_matrix <- as.matrix(Elev, wide=T)
Slope <- rast("~/data/ArcDEM/mean_slope_masked_50_res.tif")
Slope_matrix <- as.matrix(Slope, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
sl <- inla.matrix2vector(Slope_matrix)
e <- inla.matrix2vector(Elev_matrix)
table(is.na(Slope_matrix))
table(is.na(Shannon_matrix))
table(is.na(Elev_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, topo=e, slope = sl, node=node)
data <- data[-c(which(is.na(s))),]
hist(data$slope)
formula= Shannon_index ~ 1 + slope + topo +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_ele_slope_res_500 <- inla(formula,
                                        family = "gamma",
                                        data = data,
                                        control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_ele_slope_res_500)

save(Gamma_shannon_ele_slope_res_500, file="~/data/output/INLA_modelling/Gamma_shannon_ele_slope_res_500.Rdata")
get(load("~/data/output/INLA_modelling/Gamma_shannon_ele_slope_res_500.Rdata"))

## 1000m
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 100, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Elev <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 100, "_res.tif", sep=""))
Elev_matrix <- as.matrix(Elev, wide=T)
Slope <- rast("~/data/ArcDEM/mean_slope_masked_100_res.tif")
Slope_matrix <- as.matrix(Slope, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
sl <- inla.matrix2vector(Slope_matrix)
e <- inla.matrix2vector(Elev_matrix)
table(is.na(Slope_matrix))
table(is.na(Shannon_matrix))
table(is.na(Elev_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, topo=e, slope = sl, node=node)
data <- data[-c(which(is.na(s))),]
hist(data$slope)
formula= Shannon_index ~ 1 + slope + topo +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_ele_slope_res_1000 <- inla(formula,
                                        family = "gamma",
                                        data = data,
                                        control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_ele_slope_res_1000)

save(Gamma_shannon_ele_slope_res_1000, file="~/data/output/INLA_modelling/Gamma_shannon_ele_slope_res_1000.Rdata")
get(load("~/data/output/INLA_modelling/Gamma_shannon_ele_slope_res_500.Rdata"))


# plotting
Slope_df <- rbind(Gamma_shannon_ele_slope_res_100$summary.fixed[2,c(1,3,5)],
Gamma_shannon_ele_slope_res_200$summary.fixed[2,c(1,3,5)],
Gamma_shannon_ele_slope_res_300$summary.fixed[2,c(1,3,5)],
Gamma_shannon_ele_slope_res_500$summary.fixed[2,c(1,3,5)],
Gamma_shannon_ele_slope_res_1000$summary.fixed[2,c(1,3,5)])
est <- c("100m resolution", "200m resolution", "300m resolution", "500m resolution", "1000m resolution")
Slope_df <- bind_cols(as.factor(est), Slope_df)
Slope_df$model <- "model_slope"

sd_topo_df <- rbind(Gamma_shannon_ele_slope_res_100$summary.fixed[3,c(1,3,5)],
                  Gamma_shannon_ele_slope_res_200$summary.fixed[3,c(1,3,5)],
                  Gamma_shannon_ele_slope_res_300$summary.fixed[3,c(1,3,5)],
                  Gamma_shannon_ele_slope_res_500$summary.fixed[3,c(1,3,5)],
                  Gamma_shannon_ele_slope_res_1000$summary.fixed[3,c(1,3,5)])
est <- c("100m resolution", "200m resolution", "300m resolution", "500m resolution", "1000m resolution")
sd_topo_df <- bind_cols(as.factor(est), sd_topo_df)
sd_topo_df$model <- "model_topo"
slope_n_topo <- rbind(sd_topo_df, Slope_df)
colnames(slope_n_topo)[3:4] <- c("lower", "upper")

ggplot(data = slope_n_topo, 
       aes(x = mean, y = ...1, xmin = lower, xmax = upper, 
           color = model)) +
  geom_pointrange() +
  labs(title = "Model Estimates of Brain and Body Weight on REM Sleep",
       x = "Coefficient Estimate",
       y = "Predictor",
       caption = "Models fit with OLS. Error bars show the 95% confidence interval.") 
# +
#   scale_y_discrete(labels = c("Intercept", "Body Weight", "Brain Weight"))

sd_topo_plot <- slope_n_topo %>% 
  mutate(est= factor(est, levels= c("100m resolution", "200m resolution", "300m resolution", "500m resolution", "1000m resolution"))) %>% 
  ggplot() + 
  geom_pointrange(aes(x = factor(est), y = mean, ymin = lower, ymax = upper), position = position_dodge(0.5)) +
  xlab("") +
  ylab("")+
  coord_flip()+
  geom_vline(xintercept = 0, linetype='longdash', col="red")+
  theme_bw()

sd_topo_plot

sd_topo_plot_trans <- Gamma_shannon_topo_df_trans %>% 
  mutate(est= factor(est, levels= c("100m resolution", "200m resolution", "300m resolution", "500m resolution", "1000m resolution"))) %>% 
  ggplot() + 
  geom_pointrange(aes(x = factor(est), y = mean, ymin = lower, ymax = upper), position = position_dodge(0.5)) +
  xlab("coefficient sd(topography)") +
  ylab("Estimate")+
  coord_flip()+
  theme_bw()+
  geom_hline(yintercept = as.numeric(1), linetype='longdash', col="red", alpha=0.5)
sd_topo_plot_trans
geom_vline(data=ev,aes(xintercept=as.numeric(dt)))

par(mfrow=c(1,2))
sd_topo_plot
sd_topo_plot_trans
plot_grid(sd_topo_plot, sd_topo_plot_trans)


################
# model with max elevation
################
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 10, "PC1278", sep="_"))
max_ele <- rast("~/data/ArcDEM/max_ele_masked_10_res.tif")
Shannon_matrix <- as.matrix(Shannon, wide=T)
max_ele_matrix <- as.matrix(max_ele, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
max <- inla.matrix2vector(max_ele_matrix)
table(is.na(max_ele_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, max = max, node=node)
data <- data[-c(which(is.na(s))),]

formula= Shannon_index ~ 1 + max +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)

Gamma_shannon_max_ele_res_100 <- inla(formula,     
                                    family = "gamma",
                                    data = data,
                                    control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_max_ele_res_100)

glm_res100_shannon_maxele_identity <- glm(Shannon_index ~ 1 + max_ele, family = Gamma(link = "identity"), data = data)
summary(glm_res100_shannon_maxele_identity)

summary(glm(Shannon_index~1+topo+slope, family = Gamma(link = "log"), data = data))


################
# model with mean slope only
################
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 10, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast("~/data/ArcDEM/mean_slope_masked_10_res.tif")
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
save(Gamma_shannon_slope_res_100, file="~/scratch/INLA_modelling/Gamma_shannon_mean_slope_res_100.Rdata")
get(load("~/scratch/INLA_modelling/Gamma_shannon_mean_slope_res_100.Rdata"))

observed <- data$Shannon_index
fit <- Gamma_shannon_sd_slope_res_100$summary.fitted.values$mean[1:length(observed)]
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
  # coord_cartesian(xlim = c(1, 3.5))+
  scale_x_continuous(breaks=c(1,2,3), limits =c(1, 3.5))+
  scale_y_continuous(breaks=c(1,2,3), limits =c(1, 4))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
        panel.grid.major.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=margin(t = 1, r = 4, b = 1, l = 1, unit = "pt"), 
        panel.background = element_rect(fill=NA),
        text = element_text(size=20))+
  draw_text(text=paste("pseudo-R2 =", pr), x=1.7, y=3, size=20)
qq100

-sum(log(Gamma_shannon_sd_slope_res_100$cpo$cpo))

# 200m
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 20, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast("~/data/ArcDEM/mean_slope_masked_20_res.tif")
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
save(Gamma_shannon_slope_res_200, file="~/data/output/INLA_modelling/Gamma_shannon_mean_slope_res_200.Rdata")
get(load("~/scratch/INLA_modelling/Gamma_shannon_mean_slope_res_200.Rdata"))

observed <- data$Shannon_index
fit <- Gamma_shannon_slope_res_200$summary.fitted.values$mean[1:length(observed)]
# plot_inla_residuals(Gamma_shannon_slope_res_200, observed=observed)
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
  # coord_cartesian(xlim = c(1, 3.5))+
  scale_x_continuous(breaks=c(1,2,3), limits =c(1, 3.5))+
  scale_y_continuous(breaks=c(1,2,3), limits =c(1, 4))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
        panel.grid.major.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=margin(t = 1, r = 4, b = 1, l = 1, unit = "pt"), 
        panel.background = element_rect(fill=NA),
        text = element_text(size=20))+
  draw_text(text=paste("pseudo-R2 =", pr), x=1.7, y=3, size=20)
qq200
-sum(log(Gamma_shannon_slope_res_200$cpo$cpo))

# 300m
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 30, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast("~/data/ArcDEM/mean_slope_masked_30_res.tif")
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
save(Gamma_shannon_slope_res_300, file="~/data/output/INLA_modelling/Gamma_shannon_mean_slope_res_300.Rdata")
get(load("~/scratch/INLA_modelling/Gamma_shannon_mean_slope_res_300.Rdata"))

observed <- data$Shannon_index
fit <- Gamma_shannon_slope_res_300$summary.fitted.values$mean[1:length(observed)]
# plot_inla_residuals(Gamma_shannon_slope_res_200, observed=observed)
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
  # coord_cartesian(xlim = c(1, 3.5))+
  scale_x_continuous(breaks=c(1,2,3), limits =c(1, 3.5))+
  scale_y_continuous(breaks=c(1,2,3), limits =c(1, 4))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
        panel.grid.major.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=margin(t = 1, r = 4, b = 1, l = 1, unit = "pt"), 
        panel.background = element_rect(fill=NA),
        text = element_text(size=20))+
  draw_text(text=paste("pseudo-R2 =", pr), x=1.7, y=3, size=20)
qq300

-sum(log(Gamma_shannon_slope_res_300$cpo$cpo))

# 500
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 50, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast("~/data/ArcDEM/mean_slope_masked_50_res.tif")
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
save(Gamma_shannon_slope_res_500, file="~/data/output/INLA_modelling/Gamma_shannon_mean_slope_res_500.Rdata")
get(load("~/scratch/INLA_modelling/Gamma_shannon_mean_slope_res_500.Rdata"))

observed <- data$Shannon_index
fit <- Gamma_shannon_slope_res_500$summary.fitted.values$mean[1:length(observed)]
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
qq500 <- ggplot()+  
  geom_abline(slope=1, intercept=0, col="grey")+
  geom_point(data=data.frame(observed=observed, fit=fit), aes(fit, observed), pch=21)+
  ylab("")+
  xlab("")+
  # coord_cartesian(xlim = c(1, 3.5))+
  scale_x_continuous(breaks=c(1,2,3), limits =c(1, 3.5))+
  scale_y_continuous(breaks=c(1,2,3), limits =c(1, 4))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
        panel.grid.major.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=margin(t = 1, r = 4, b = 1, l = 1, unit = "pt"), 
        panel.background = element_rect(fill=NA),
        text = element_text(size=20))+
  draw_text(text=paste("pseudo-R2 =", pr), x=1.7, y=3, size=20)
qq500

-sum(log(Gamma_shannon_slope_res_500$cpo$cpo))


# 1000m
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 100, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast("~/data/ArcDEM/mean_slope_masked_100_res.tif")
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
save(Gamma_shannon_slope_res_1000, file="~/data/output/INLA_modelling/Gamma_shannon_mean_slope_res_1000.Rdata")
get(load("~/scratch/INLA_modelling/Gamma_shannon_mean_slope_res_1000.Rdata"))

observed <- data$Shannon_index
fit <- Gamma_shannon_slope_res_1000$summary.fitted.values$mean[1:length(observed)]
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
qq1000 <- ggplot()+  
  geom_abline(slope=1, intercept=0, col="grey")+
  geom_point(data=data.frame(observed=observed, fit=fit), aes(fit, observed), pch=21)+
  ylab("")+
  xlab("")+
  # coord_cartesian(xlim = c(1, 3.5))+
  scale_x_continuous(breaks=c(1,2,3), limits =c(1, 3.6))+
  scale_y_continuous(breaks=c(1,2,3), limits =c(1, 4))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
        panel.grid.major.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=margin(t = 1, r = 4, b = 1, l = 1, unit = "pt"), 
        panel.background = element_rect(fill=NA),
        # axis.title.x = element_text(size = 20)
        text = element_text(size=20))+
  draw_text(text=paste("pseudo-R2 =", pr), x=1.7, y=3, size=20)
qq1000

-sum(log(Gamma_shannon_slope_res_1000$cpo$cpo))

qqplot_sd_slope_together <- plot_grid(qq100, qq200, qq300, qq500, qq1000, align = "v", ncol=1, axis = "r", labels="AUTO")
qqplot_sd_slope_together
save_plot(qqplot_sd_slope_together, filename = "~/data/output/final_plot/qqplot_sd_slope_together.png", base_height = 25, base_width = 10,bg="white")

library(grid)
library(gridExtra)
x.grob <- textGrob(paste("Fitted Shannon index from model with", "standard deviation of topographic slope as predictor", sep="\n"), gp=gpar(fontsize=25))
y.grob <- textGrob("Observe Shannon index", rot=90, gp=gpar(fontsize=25)) 
qqplot_sd_slope_together_final <- grid.arrange(arrangeGrob(qqplot_sd_slope_together,bottom = x.grob, left = y.grob))
qqplot_sd_slope_together_final
save_plot(qqplot_sd_slope_together_final, filename = "~/data/output/final_plot/qqplot_sd_slope_together_final.png", base_height = 25, base_width = 10,bg="white")

# plotting

Slope_df <- rbind(Gamma_shannon_slope_res_100$summary.fixed[2,c(1,3,5)],
                  Gamma_shannon_slope_res_200$summary.fixed[2,c(1,3,5)],
                  Gamma_shannon_slope_res_300$summary.fixed[2,c(1,3,5)],
                  Gamma_shannon_slope_res_500$summary.fixed[2,c(1,3,5)],
                  Gamma_shannon_slope_res_1000$summary.fixed[2,c(1,3,5)])
est <- c("100m resolution", "200m resolution", "300m resolution", "500m resolution", "1000m resolution")
Slope_df <- cbind(as.factor(est), Slope_df)
colnames(Slope_df)[c(1,3:4)] <- c("resolution","lower", "upper")
Slope_df[,2:4] <- exp(Slope_df[,2:4])
Slope_df$resolution <- factor(Slope_df$resolution, levels = c("100m resolution", "200m resolution", "300m resolution", "500m resolution", "1000m resolution"))

mean_slope_plot_trans <- ggplot(data = Slope_df, 
       aes(x = mean, y = resolution, xmin = lower, xmax = upper)) +
  geom_pointrange() +
  labs(title = "Model Estimates of mean slope on estimated Shannon index",
       x = "Coefficient Estimate",
       y = "resolution",
       caption = "Models fit with INLA. Error bars show the 95% confidence interval.")+
  coord_cartesian(xlim = c(0.982, 1))+
  geom_vline(xintercept = as.numeric(1), linetype='longdash', col="red", alpha=0.5)+
  theme_bw()
  


#################
# Model with coefficient of variation
#################
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 10, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Elev <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 10, "_res.tif", sep=""))
Elev_matrix <- as.matrix(Elev, wide=T)
Mean_ele <- rast(paste("~/data/ArcDEM/mean_ele_masked_", 10, "_res.tif", sep=""))
Mean_ele_matrix <- as.matrix(Mean_ele, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
e <- inla.matrix2vector(Elev_matrix)
me <- inla.matrix2vector(Mean_ele_matrix)
table(is.na(Shannon_matrix))
table(is.na(Elev_matrix))
table(is.na(Mean_ele_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, topo=e, mean_ele = me, node=node)
data <- data[-c(which(is.na(s))),]
data$coeff_var <- data$topo/data$mean_ele
# data <- data[-11798,] 
plot(data$Shannon_index~data$coeff_var)

# formula= Shannon_index ~ 1 + coeff_var +
#   f(node, model="matern2d", nrow=nrow, ncol=ncol)
# Gamma_shannon_coeffvar_res_100 <- inla(formula,     
#                                      family = "gamma",
#                                      data = data,
#                                      control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE, config=T), verbose=TRUE)
# 
# summary(Gamma_shannon_coeffvar_res_100)
# save(Gamma_shannon_coeffvar_res_100, file="~/scratch/INLA_modelling/Gamma_shannon_coeffvar_res_100.Rdata")
get(load("~/scratch/INLA_modelling/Gamma_shannon_coeffvar_res_100.Rdata"))
# coeff_var_glm_100 <- glm(Shannon_index ~ coeff_var, family = Gamma(link = "identity"), data = data)

observed <- data$Shannon_index
fit <- Gamma_shannon_coeffvar_res_100$summary.fitted.values$mean[1:length(observed)]
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
  # coord_cartesian(xlim = c(1, 3.5))+
  scale_x_continuous(breaks=c(1,2,3), limits =c(1, 3.5))+
  scale_y_continuous(breaks=c(1,2,3), limits =c(1, 4))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
        panel.grid.major.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=margin(t = 1, r = 4, b = 1, l = 1, unit = "pt"), 
        panel.background = element_rect(fill=NA),
        text = element_text(size=20))+
  draw_text(text=paste("pseudo-R2 =", pr), x=1.7, y=3, size=20)
qq100
-sum(log(Gamma_shannon_coeffvar_res_100$cpo$cpo))


# formula= Shannon_index ~ 1 + mean_ele +
#   f(node, model="matern2d", nrow=nrow, ncol=ncol)
# Gamma_shannon_mean_ele_res_100 <- inla(formula,     
#                                        family = "gamma",
#                                        data = data,
#                                        control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
# save(Gamma_shannon_mean_ele_res_100, file="~/data/output/INLA_modelling/Gamma_shannon_mean_ele_res_100.Rdata")
# get(load("~/data/output/INLA_modelling/Gamma_shannon_mean_ele_res_100.Rdata"))



# 200
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 20, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Elev <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 20, "_res.tif", sep=""))
Elev_matrix <- as.matrix(Elev, wide=T)
Mean_ele <- Elev <- rast(paste("~/data/ArcDEM/mean_ele_masked_", 20, "_res.tif", sep=""))
Mean_ele_matrix <- as.matrix(Mean_ele, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
e <- inla.matrix2vector(Elev_matrix)
me <- inla.matrix2vector(Mean_ele_matrix)
table(is.na(Shannon_matrix))
table(is.na(Elev_matrix))
table(is.na(Mean_ele_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, topo=e, mean_ele = me, node=node)
data <- data[-c(which(is.na(s))),]
data$coeff_var <- data$topo/data$mean_ele

# formula= Shannon_index ~ 1 + coeff_var +
#   f(node, model="matern2d", nrow=nrow, ncol=ncol)
# Gamma_shannon_coeffvar_res_200 <- inla(formula,     
#                                        family = "gamma",
#                                        data = data,
#                                        control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE, config=T), verbose=TRUE)
# 
# summary(Gamma_shannon_coeffvar_res_200)
# save(Gamma_shannon_coeffvar_res_200, file="~/scratch/INLA_modelling/Gamma_shannon_coeffvar_res_200.Rdata")
get(load("~/scratch/INLA_modelling/Gamma_shannon_coeffvar_res_200.Rdata"))

observed <- data$Shannon_index
fit <- Gamma_shannon_coeffvar_res_200$summary.fitted.values$mean[1:length(observed)]
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
  # coord_cartesian(xlim = c(1, 3.5))+
  scale_x_continuous(breaks=c(1,2,3), limits =c(1, 3.5))+
  scale_y_continuous(breaks=c(1,2,3), limits =c(1, 4))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
        panel.grid.major.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=margin(t = 1, r = 4, b = 1, l = 1, unit = "pt"), 
        panel.background = element_rect(fill=NA),
        text = element_text(size=20))+
  draw_text(text=paste("pseudo-R2 =", pr), x=1.7, y=3, size=20)
qq200
-sum(log(Gamma_shannon_coeffvar_res_200$cpo$cpo))

# 
# formula= Shannon_index ~ 1 + mean_ele +
#   f(node, model="matern2d", nrow=nrow, ncol=ncol)
# Gamma_shannon_mean_ele_res_200 <- inla(formula,     
#                                        family = "gamma",
#                                        data = data,
#                                        control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
# save(Gamma_shannon_mean_ele_res_200, file="~/data/output/INLA_modelling/Gamma_shannon_mean_ele_res_200.Rdata")
# get(load("~/data/output/INLA_modelling/Gamma_shannon_mean_ele_res_200.Rdata"))


# 300
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 30, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Elev <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 30, "_res.tif", sep=""))
Elev_matrix <- as.matrix(Elev, wide=T)
Mean_ele <- Elev <- rast(paste("~/data/ArcDEM/mean_ele_masked_", 30, "_res.tif", sep=""))
Mean_ele_matrix <- as.matrix(Mean_ele, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
e <- inla.matrix2vector(Elev_matrix)
me <- inla.matrix2vector(Mean_ele_matrix)
table(is.na(Shannon_matrix))
table(is.na(Elev_matrix))
table(is.na(Mean_ele_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, topo=e, mean_ele = me, node=node)
data <- data[-c(which(is.na(s))),]
data$coeff_var <- data$topo/data$mean_ele

formula= Shannon_index ~ 1 + coeff_var +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_coeffvar_res_300 <- inla(formula,
                                       family = "gamma",
                                       data = data,
                                       control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE, config=T), verbose=TRUE)

# summary(Gamma_shannon_coeffvar_res_300)
# save(Gamma_shannon_coeffvar_res_300, file="~/scratch/INLA_modelling/Gamma_shannon_coeffvar_res_300.Rdata")
get(load("~/scratch/INLA_modelling/Gamma_shannon_coeffvar_res_300.Rdata"))

observed <- data$Shannon_index
fit <- Gamma_shannon_coeffvar_res_300$summary.fitted.values$mean[1:length(observed)]
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
  # coord_cartesian(xlim = c(1, 3.5))+
  scale_x_continuous(breaks=c(1,2,3), limits =c(1, 3.5))+
  scale_y_continuous(breaks=c(1,2,3), limits =c(1, 4))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
        panel.grid.major.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=margin(t = 1, r = 4, b = 1, l = 1, unit = "pt"), 
        panel.background = element_rect(fill=NA),
        text = element_text(size=20))+
  draw_text(text=paste("pseudo-R2 =", pr), x=1.7, y=3, size=20)
qq300
-sum(log(Gamma_shannon_coeffvar_res_300$cpo$cpo))

# 
# formula= Shannon_index ~ 1 + mean_ele +
#   f(node, model="matern2d", nrow=nrow, ncol=ncol)
# Gamma_shannon_mean_ele_res_300 <- inla(formula,     
#                                        family = "gamma",
#                                        data = data,
#                                        control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
# save(Gamma_shannon_mean_ele_res_300, file="~/data/output/INLA_modelling/Gamma_shannon_mean_ele_res_300.Rdata")
# get(load("~/data/output/INLA_modelling/Gamma_shannon_mean_ele_res_300.Rdata"))


# 500
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 50, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Elev <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 50, "_res.tif", sep=""))
Elev_matrix <- as.matrix(Elev, wide=T)
Mean_ele <- Elev <- rast(paste("~/data/ArcDEM/mean_ele_masked_", 50, "_res.tif", sep=""))
Mean_ele_matrix <- as.matrix(Mean_ele, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
e <- inla.matrix2vector(Elev_matrix)
me <- inla.matrix2vector(Mean_ele_matrix)
table(is.na(Shannon_matrix))
table(is.na(Elev_matrix))
table(is.na(Mean_ele_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, topo=e, mean_ele = me, node=node)
data <- data[-c(which(is.na(s))),]
data$coeff_var <- data$topo/data$mean_ele

# formula= Shannon_index ~ 1 + coeff_var +
#   f(node, model="matern2d", nrow=nrow, ncol=ncol)
# Gamma_shannon_coeffvar_res_500 <- inla(formula,     
#                                        family = "gamma",
#                                        data = data,
#                                        control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE, config=T), verbose=TRUE)
# 
# summary(Gamma_shannon_coeffvar_res_500)
# save(Gamma_shannon_coeffvar_res_500, file="~/scratch/INLA_modelling/Gamma_shannon_coeffvar_res_500.Rdata")
get(load("~/scratch/INLA_modelling/Gamma_shannon_coeffvar_res_500.Rdata"))

observed <- data$Shannon_index
fit <- Gamma_shannon_coeffvar_res_500$summary.fitted.values$mean[1:length(observed)]
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
qq500 <- ggplot()+  
  geom_abline(slope=1, intercept=0, col="grey")+
  geom_point(data=data.frame(observed=observed, fit=fit), aes(fit, observed), pch=21)+
  ylab("")+
  xlab("")+
  # coord_cartesian(xlim = c(1, 3.5))+
  scale_x_continuous(breaks=c(1,2,3), limits =c(1, 3.5))+
  scale_y_continuous(breaks=c(1,2,3), limits =c(1, 4))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
        panel.grid.major.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=margin(t = 1, r = 4, b = 1, l = 1, unit = "pt"), 
        panel.background = element_rect(fill=NA),
        text = element_text(size=20))+
  draw_text(text=paste("pseudo-R2 =", pr), x=1.7, y=3, size=20)
qq500
-sum(log(Gamma_shannon_coeffvar_res_500$cpo$cpo))

# formula= Shannon_index ~ 1 + mean_ele +
#   f(node, model="matern2d", nrow=nrow, ncol=ncol)
# Gamma_shannon_mean_ele_res_500 <- inla(formula,     
#                                        family = "gamma",
#                                        data = data,
#                                        control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
# save(Gamma_shannon_mean_ele_res_500, file="~/data/output/INLA_modelling/Gamma_shannon_mean_ele_res_500.Rdata")
# get(load("~/data/output/INLA_modelling/Gamma_shannon_mean_ele_res_500.Rdata"))


# 1000
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 100, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Elev <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 100, "_res.tif", sep=""))
Elev_matrix <- as.matrix(Elev, wide=T)
Mean_ele <- Elev <- rast(paste("~/data/ArcDEM/mean_ele_masked_", 100, "_res.tif", sep=""))
Mean_ele_matrix <- as.matrix(Mean_ele, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
e <- inla.matrix2vector(Elev_matrix)
me <- inla.matrix2vector(Mean_ele_matrix)
table(is.na(Shannon_matrix))
table(is.na(Elev_matrix))
table(is.na(Mean_ele_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, topo=e, mean_ele = me, node=node)
data <- data[-c(which(is.na(s))),]
data$coeff_var <- data$topo/data$mean_ele

formula= Shannon_index ~ 1 + coeff_var +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_coeffvar_res_1000 <- inla(formula,
                                       family = "gamma",
                                       data = data,
                                       control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE, config=T), verbose=TRUE)

# summary(Gamma_shannon_coeffvar_res_1000)
# save(Gamma_shannon_coeffvar_res_1000, file="~/scratch/INLA_modelling/Gamma_shannon_coeffvar_res_1000.Rdata")
get(load("~/scratch/INLA_modelling/Gamma_shannon_coeffvar_res_1000.Rdata"))

observed <- data$Shannon_index
fit <- Gamma_shannon_coeffvar_res_1000$summary.fitted.values$mean[1:length(observed)]
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
qq1000 <- ggplot()+  
  geom_abline(slope=1, intercept=0, col="grey")+
  geom_point(data=data.frame(observed=observed, fit=fit), aes(fit, observed), pch=21)+
  ylab("")+
  xlab("")+
  # coord_cartesian(xlim = c(1, 3.5))+
  scale_x_continuous(breaks=c(1,2,3), limits =c(1, 3.5))+
  scale_y_continuous(breaks=c(1,2,3), limits =c(1, 4))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
        panel.grid.major.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=margin(t = 1, r = 4, b = 1, l = 1, unit = "pt"), 
        panel.background = element_rect(fill=NA),
        text = element_text(size=20))+
  draw_text(text=paste("pseudo-R2 =", pr), x=1.7, y=3, size=20)
qq1000
-sum(log(Gamma_shannon_coeffvar_res_1000$cpo$cpo))


# formula= Shannon_index ~ 1 + mean_ele +
#   f(node, model="matern2d", nrow=nrow, ncol=ncol)
# Gamma_shannon_mean_ele_res_1000 <- inla(formula,     
#                                        family = "gamma",
#                                        data = data,
#                                        control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
# save(Gamma_shannon_mean_ele_res_1000, file="~/data/output/INLA_modelling/Gamma_shannon_mean_ele_res_1000.Rdata")
# get(load("~/data/output/INLA_modelling/Gamma_shannon_mean_ele_res_1000.Rdata"))

qqplot_CV_ele_together <- plot_grid(qq100, qq200, qq300, qq500, qq1000, align = "v", ncol=1, axis = "r", labels="AUTO", label_size = 25)
qqplot_CV_ele_together
save_plot(qqplot_CV_ele_together, filename = "~/data/output/final_plot/qqplot_CV_ele_together.png", base_height = 25, base_width = 10,bg="white")

library(grid)
library(gridExtra)
x.grob <- textGrob(paste("Fitted Shannon index from model with", "coefficient of variation of elevation as predictor", sep="\n"), gp=gpar(fontsize=25), vjust=0.2)
y.grob <- textGrob("Observe Shannon index", rot=90, gp=gpar(fontsize=25)) 
qqplot_CV_ele_together_final <- grid.arrange(arrangeGrob(qqplot_CV_ele_together,bottom = x.grob, left = y.grob))
qqplot_CV_ele_together_final
save_plot(qqplot_CV_ele_together_final, filename = "~/data/output/final_plot/qqplot_CV_ele_together_final.png", base_height = 25, base_width = 10,bg="white")

# plotting CV results
CV_df <- rbind(Gamma_shannon_coeffvar_res_100$summary.fixed[2,c(1,3,5)],
               Gamma_shannon_coeffvar_res_200$summary.fixed[2,c(1,3,5)],
               Gamma_shannon_coeffvar_res_300$summary.fixed[2,c(1,3,5)],
               Gamma_shannon_coeffvar_res_500$summary.fixed[2,c(1,3,5)],
               Gamma_shannon_coeffvar_res_1000$summary.fixed[2,c(1,3,5)])
est <- c("100m resolution", "200m resolution", "300m resolution", "500m resolution", "1000m resolution")
CV_df <- cbind(as.factor(est), CV_df)
colnames(CV_df)[c(1,3:4)] <- c("resolution","lower", "upper")
# Slope_df[,2:4] <- exp(Slope_df[,2:4])
CV_df$resolution <- factor(CV_df$resolution, levels = c("100m resolution", "200m resolution", "300m resolution", "500m resolution", "1000m resolution"))

CV_coeff_plot <- ggplot(data = CV_df, 
                                aes(x = mean, y = resolution, xmin = lower, xmax = upper)) +
  geom_pointrange(fatten=2.5) +
  labs(
    # title = "Model Estimates of the coefficient of variation in elevation on Shannon index estimates",
       x = "Coefficient Estimate",
       y = "resolution"
       # ,caption = "Models fit with INLA using a log-link function. Error bars show the 95% confidence interval."
       )+
  geom_vline(xintercept = as.numeric(0), linetype='longdash', col="red", alpha=0.5)+
  coord_cartesian(xlim = c(-0.01, 0.01))+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
CV_coeff_plot


plot(data$topo~data$mean_ele)
mean_ele_sd_ele <- ggplot(data, aes(x=mean_ele, y=topo))+
  geom_point(shape=21)+
  theme_cowplot()+
  labs(x="Mean elevation", y=expression(paste(sigma, " (elevation)")))+
  theme(axis.line=element_line(linewidth =0.5), axis.text.y = element_text(size=12), axis.text.x = element_text(size=12), axis.title.x = element_text(size=14), 
        axis.title.y=element_text(size=14), panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
mean_ele_sd_ele

save_plot(mean_ele_sd_ele, filename = "~/data/output/final_plot/scatterplot_mean_ele_sd_ele.png", bg="white", base_height = 4, base_asp = 1.4)

#######
# link topography and sd(slope)
#######
topo_10 <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 10, "_res.tif", sep=""))
slope_10 <- rast("~/data/ArcDEM/mean_slope_masked_10_res.tif")
topo_20 <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 20, "_res.tif", sep=""))
slope_20 <- rast("~/data/ArcDEM/mean_slope_masked_20_res.tif")
topo_30 <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 30, "_res.tif", sep=""))
slope_30 <- rast("~/data/ArcDEM/mean_slope_masked_30_res.tif")
topo_50 <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 50, "_res.tif", sep=""))
slope_50 <- rast("~/data/ArcDEM/mean_slope_masked_50_res.tif")
topo_100 <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 100, "_res.tif", sep=""))
slope_100 <- rast("~/data/ArcDEM/mean_slope_masked_100_res.tif")
plot(topo_10)
plot(slope_10)
plot(topo_20)
plot(slope_20)
plot(topo_30)
plot(slope_30)
plot(topo_50)
plot(slope_50)
plot(topo_100)
plot(slope_100)



##############
# Model with sd(slope), including glm plotting
##############
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 10, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast("~/data/ArcDEM/sd_slope_masked_10_res.tif")
Slope_matrix <- as.matrix(Slope, wide=T)
mean_slope <- rast("~/data/ArcDEM/mean_slope_masked_10_res.tif")
mean_slope_matrix <- as.matrix(mean_slope)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
sl <- inla.matrix2vector(Slope_matrix)
msl <- inla.matrix2vector(mean_slope_matrix)
table(is.na(Slope_matrix))
table(is.na(Shannon_matrix))
table(is.na(mean_slope_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, slope = sl, mean_sl = msl, node=node)
data <- data[-c(which(is.na(s))),]
data$slope_CV <- data$slope/data$mean_sl

# formula= Shannon_index ~ 1 + slope +
#   f(node, model="matern2d", nrow=nrow, ncol=ncol)
# Gamma_shannon_sd_slope_res_100 <- inla(formula,     
#                                     family = "gamma",
#                                     data = data,
#                                     control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE, config=T), 
#                                     control.predictor=list(compute=T, link=1), 
#                                     verbose=TRUE)
# summary(Gamma_shannon_sd_slope_res_100)
# save(Gamma_shannon_sd_slope_res_100, file="~/scratch/INLA_modelling/Gamma_shannon_sd_slope_res_100.Rdata")
get(load("~/scratch/INLA_modelling/Gamma_shannon_sd_slope_res_100.Rdata"))

observed <- data$Shannon_index
fit <- Gamma_shannon_sd_slope_res_100$summary.fitted.values$mean[1:length(observed)]
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
  # coord_cartesian(xlim = c(1, 3.5))+
  scale_x_continuous(breaks=c(1,2,3), limits =c(1, 3.5))+
  scale_y_continuous(breaks=c(1,2,3), limits =c(1, 4))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
        panel.grid.major.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=margin(t = 1, r = 4, b = 1, l = 1, unit = "pt"), 
        panel.background = element_rect(fill=NA),
        text = element_text(size=20))+
  draw_text(text=paste("pseudo-R2 =", pr), x=1.7, y=3, size=20)
qq100

-sum(log(Gamma_shannon_sd_slope_res_100$cpo$cpo))

# # glm plotting
# res_100_glm <- glm(Shannon_index ~ 1+slope, family = Gamma(link="log"), data = data)
# summary(res_100_glm)
# plot(res_100_glm)
# res <- residuals(res_100_glm, type="partial")
# plot(data$slope, res)
# plot(fitted(res_100_glm), residuals(res_100_glm, type="pearson"))
# (res_100_plot <- ggplot(
#   data,
#   aes(x = slope, y = Shannon_index)) +
#     geom_point(shape=21) +
#     scale_y_continuous(limits = c(1, 4)) +
#     labs(x = "sd(slope)", y="Shannon index", title = "100m resolution")+
#     theme_cowplot())
# 
# res_100_pred <- predict(
#   res_100_glm,
#   newdata = data.frame(slope = seq(0, 15, 0.01)),
#   type = "response"
# )
# 
# res_100_pred_df <- data.frame(
#   slope = seq(0,15,0.01),
#   Shannon_index = res_100_pred
# )
# res_100_plot_pred <- res_100_plot +
#   geom_line(data = res_100_pred_df, 
#             colour = "blue")
# res_100_plot_pred
# 
# # with coeff variation
# formula= Shannon_index ~ 1 + slope_CV +
#   f(node, model="matern2d", nrow=nrow, ncol=ncol)
# Gamma_shannon_slopeCV_res_100 <- inla(formula,     
#                                        family = "gamma",
#                                        data = data,
#                                        control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
# summary(Gamma_shannon_slopeCV_res_100)

# corr <- raster.modified.ttest(Shannon, Slope)


# 200m
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 20, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast("~/data/ArcDEM/sd_slope_masked_20_res.tif")
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

# formula= Shannon_index ~ 1 + slope +
#   f(node, model="matern2d", nrow=nrow, ncol=ncol)
# Gamma_shannon_sd_slope_res_200 <- inla(formula,     
#                                     family = "gamma",
#                                     data = data,
#                                     control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE, config=T), verbose=TRUE)
# summary(Gamma_shannon_sd_slope_res_200)
# save(Gamma_shannon_sd_slope_res_200, file="~/data/output/INLA_modelling/Gamma_shannon_sd_slope_res_200.Rdata")
get(load("~/scratch/INLA_modelling/Gamma_shannon_sd_slope_res_200.Rdata"))

observed <- data$Shannon_index
fit <- Gamma_shannon_sd_slope_res_200$summary.fitted.values$mean[1:length(observed)]
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
  # coord_cartesian(xlim = c(1, 3.5))+
  scale_x_continuous(breaks=c(1,2,3), limits =c(1, 3.5))+
  scale_y_continuous(breaks=c(1,2,3), limits =c(1, 4))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
        panel.grid.major.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=margin(t = 1, r = 4, b = 1, l = 1, unit = "pt"), 
        panel.background = element_rect(fill=NA),
        text = element_text(size=20))+
  draw_text(text=paste("pseudo-R2 =", pr), x=1.7, y=3, size=20)
qq200

-sum(log(Gamma_shannon_sd_slope_res_200$cpo$cpo))

# # glm plotting
# res_200_glm <- glm(Shannon_index ~ slope, family = Gamma(link = "log"), data = data)
# summary(res_200_glm)
# (res_200_plot <- ggplot(
#   data,
#   aes(x = slope, y = Shannon_index)) +
#     geom_point(shape=21) +
#     scale_y_continuous(limits = c(1, 4)) +
#     labs(x = "sd(slope)", y="Shannon index", title = "200m resolution")+
#     theme_cowplot())
# range(data$slope)
# res_200_pred <- predict(
#   res_200_glm,
#   newdata = data.frame(slope = seq(0, 15, 0.01)),
#   type = "response"
# )
# 
# res_200_pred_df <- data.frame(
#   slope = seq(0,15,0.01),
#   Shannon_index = res_200_pred
# )
# res_200_plot_pred <- res_200_plot +
#   geom_line(data = res_200_pred_df, 
#             colour = "blue")
# res_200_plot_pred



# 300m
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 30, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast("~/data/ArcDEM/sd_slope_masked_30_res.tif")
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

# formula= Shannon_index ~ 1 + slope +
#   f(node, model="matern2d", nrow=nrow, ncol=ncol)
# Gamma_shannon_sd_slope_res_300 <- inla(formula,     
#                                     family = "gamma",
#                                     data = data,
#                                     control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE, config=T), verbose=TRUE)
# formula= Shannon_index ~ 1 + log(slope) +
#   f(node, model="matern2d", nrow=nrow, ncol=ncol)
# sn_shannon_sd_slope_res_300 <- inla(formula,     
#                                        family = "sn",
#                                        data = data,
#                                        control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE, config=T), verbose=TRUE)
# summary(sn_shannon_sd_slope_res_300)
# observed = data$Shannon_index
# plot_inla_residuals(sn_shannon_sd_slope_res_300, observed = observed)

# summary(Gamma_shannon_sd_slope_res_300)
# save(Gamma_shannon_sd_slope_res_300, file="~/data/output/INLA_modelling/Gamma_shannon_sd_slope_res_300.Rdata")
get(load("~/scratch/INLA_modelling/Gamma_shannon_sd_slope_res_300.Rdata"))

observed <- data$Shannon_index
fit <- Gamma_shannon_sd_slope_res_300$summary.fitted.values$mean[1:length(observed)]
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
  # coord_cartesian(xlim = c(1, 3.5))+
  scale_x_continuous(breaks=c(1,2,3), limits =c(1, 3.5))+
  scale_y_continuous(breaks=c(1,2,3), limits =c(1, 4))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
        panel.grid.major.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=margin(t = 1, r = 4, b = 1, l = 1, unit = "pt"), 
        panel.background = element_rect(fill=NA),
        text = element_text(size=20))+
  draw_text(text=paste("pseudo-R2 =", pr), x=1.7, y=3, size=20)
qq300
-sum(log(Gamma_shannon_sd_slope_res_300$cpo$cpo))

 
# # glm plotting
# res_300_glm <- glm(Shannon_index ~ slope, family = Gamma(link = "log"), data = data)
# summary(res_300_glm)
# (res_300_plot <- ggplot(
#   data,
#   aes(x = slope, y = Shannon_index)) +
#     geom_point(shape=21) +
#     scale_y_continuous(limits = c(1, 4)) +
#     labs(x = "sd(slope)", y="Shannon index", title = "300m resolution")+
#     theme_cowplot())
# range(data$slope)
# res_300_pred <- predict(
#   res_300_glm,
#   newdata = data.frame(slope = seq(0, 15, 0.01)),
#   type = "response"
# )
# 
# res_300_pred_df <- data.frame(
#   slope = seq(0,15,0.01),
#   Shannon_index = res_300_pred
# )
# res_300_plot_pred <- res_300_plot +
#   geom_line(data = res_300_pred_df, 
#             colour = "blue")
# res_300_plot_pred

# 500
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 50, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast("~/data/ArcDEM/sd_slope_masked_50_res.tif")
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

# formula= Shannon_index ~ 1 + slope +
#   f(node, model="matern2d", nrow=nrow, ncol=ncol)
# Gamma_shannon_sd_slope_res_500 <- inla(formula,     
#                                     family = "gamma",
#                                     data = data,
#                                     control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE, config=T), verbose=TRUE)
# summary(Gamma_shannon_sd_slope_res_500)
# save(Gamma_shannon_sd_slope_res_500, file="~/data/output/INLA_modelling/Gamma_shannon_sd_slope_res_500.Rdata")
get(load("~/scratch/INLA_modelling/Gamma_shannon_sd_slope_res_500.Rdata"))

observed <- data$Shannon_index
fit <- Gamma_shannon_sd_slope_res_500$summary.fitted.values$mean[1:length(observed)]
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
qq500 <- ggplot()+  
  geom_abline(slope=1, intercept=0, col="grey")+
  geom_point(data=data.frame(observed=observed, fit=fit), aes(fit, observed), pch=21)+
  ylab("")+
  xlab("")+
  # coord_cartesian(xlim = c(1, 3.5))+
  scale_x_continuous(breaks=c(1,2,3), limits =c(1, 3.5))+
  scale_y_continuous(breaks=c(1,2,3), limits =c(1, 4))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
        panel.grid.major.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=margin(t = 1, r = 4, b = 1, l = 1, unit = "pt"), 
        panel.background = element_rect(fill=NA),
        text = element_text(size=20))+
  draw_text(text=paste("pseudo-R2 =", pr), x=1.7, y=3, size=20)
qq500
-sum(log(Gamma_shannon_sd_slope_res_500$cpo$cpo))


# # glm plotting
# res_500_glm <- glm(Shannon_index ~ slope, family = Gamma(link = "log"), data = data)
# summary(res_500_glm)
# (res_500_plot <- ggplot(
#   data,
#   aes(x = slope, y = Shannon_index)) +
#     geom_point(shape=21) +
#     scale_y_continuous(limits = c(1, 4)) +
#     labs(x = "sd(slope)", y="Shannon index", title = "500m resolution")+
#     theme_cowplot())
# range(data$slope)
# res_500_pred <- predict(
#   res_500_glm,
#   newdata = data.frame(slope = seq(0, 15, 0.01)),
#   type = "response"
# )
# 
# res_500_pred_df <- data.frame(
#   slope = seq(0,15,0.01),
#   Shannon_index = res_500_pred
# )
# res_500_plot_pred <- res_500_plot +
#   geom_line(data = res_500_pred_df, 
#             colour = "blue")
# res_500_plot_pred

# 1000m
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 100, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast("~/data/ArcDEM/sd_slope_masked_100_res.tif")
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

# formula= Shannon_index ~ 1 + slope +
#   f(node, model="matern2d", nrow=nrow, ncol=ncol)
# Gamma_shannon_sd_slope_res_1000 <- inla(formula,     
#                                      family = "gamma",
#                                      data = data,
#                                      control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE, config=T), verbose=TRUE)
# summary(Gamma_shannon_sd_slope_res_1000)
# save(Gamma_shannon_sd_slope_res_1000, file="~/data/output/INLA_modelling/Gamma_shannon_sd_slope_res_1000.Rdata")
get(load("~/scratch/INLA_modelling/Gamma_shannon_sd_slope_res_1000.Rdata"))

observed <- data$Shannon_index
fit <- Gamma_shannon_sd_slope_res_1000$summary.fitted.values$mean[1:length(observed)]
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
qq1000 <- ggplot()+  
  geom_abline(slope=1, intercept=0, col="grey")+
  geom_point(data=data.frame(observed=observed, fit=fit), aes(fit, observed), pch=21)+
  ylab("")+
  xlab("")+
  # coord_cartesian(xlim = c(1, 3.5))+
  scale_x_continuous(breaks=c(1,2,3), limits =c(1, 3.5))+
  scale_y_continuous(breaks=c(1,2,3), limits =c(1, 4))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
        panel.grid.major.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=unit(c(1, 4, 1, 1),"pt"), 
        panel.background = element_rect(fill=NA),
        text = element_text(size=20))+
  draw_text(text=paste("pseudo-R2 =", pr), x=1.7, y=3, size=20)
qq1000
-sum(log(Gamma_shannon_sd_slope_res_1000$cpo$cpo))

# # glm plotting
# res_1000_glm <- glm(Shannon_index ~ slope, family = Gamma(link = "log"), data = data)
# summary(res_1000_glm)
# (res_1000_plot <- ggplot(
#   data,
#   aes(x = slope, y = Shannon_index)) +
#     geom_point(shape=21) +
#     scale_y_continuous(limits = c(1, 4)) +
#     labs(x = "sd(slope)", y="Shannon index", title = "1000m resolution")+
#     theme_cowplot())
# range(data$slope)
# res_1000_pred <- predict(
#   res_1000_glm,
#   newdata = data.frame(slope = seq(0, 15, 0.01)),
#   type = "response"
# )
# 
# res_1000_pred_df <- data.frame(
#   slope = seq(0,15,0.01),
#   Shannon_index = res_1000_pred
# )
# res_1000_plot_pred <- res_1000_plot +
#   geom_line(data = res_1000_pred_df, 
#             colour = "blue")
# res_1000_plot_pred

# qqplot plotting together 
qqplot_sd_slope_together <- plot_grid(qq100, qq200, qq300, qq500, qq1000, align = "v", ncol=1, axis = "r", labels="AUTO", label_size = 25)
qqplot_sd_slope_together
save_plot(qqplot_sd_slope_together, filename = "~/data/output/final_plot/qqplot_sd_slope_together.png", base_height = 25, base_width = 10,bg="white")

library(grid)
library(gridExtra)
x.grob <- textGrob(paste("Fitted Shannon index from model with", "standard variation of topographic slope as predictor", sep="\n"), gp=gpar(fontsize=25), vjust=0.2, hjust=NULL)
y.grob <- textGrob("Observe Shannon index", rot=90, gp=gpar(fontsize=25)) 
qqplot_sd_slope_together_final <- grid.arrange(arrangeGrob(qqplot_sd_slope_together, bottom = x.grob, left = y.grob))
qqplot_sd_slope_together_final
save_plot(qqplot_sd_slope_together_final, filename = "~/data/output/final_plot/qqplot_sd_slope_together_final.png", base_height = 25, base_width = 10,bg="white")



# inla coefficient plotting
Sd_slope_df <- rbind(Gamma_shannon_sd_slope_res_100$summary.fixed[2,c(1,3,5)],
                     Gamma_shannon_sd_slope_res_200$summary.fixed[2,c(1,3,5)],
                     Gamma_shannon_sd_slope_res_300$summary.fixed[2,c(1,3,5)],
                     Gamma_shannon_sd_slope_res_500$summary.fixed[2,c(1,3,5)],
                     Gamma_shannon_sd_slope_res_1000$summary.fixed[2,c(1,3,5)])
est <- c("100m resolution", "200m resolution", "300m resolution", "500m resolution", "1000m resolution")
Sd_slope_df <- cbind(as.factor(est), Sd_slope_df)
colnames(Sd_slope_df)[c(1,3:4)] <- c("resolution","lower", "upper")
# Sd_slope_df[2:4] <- exp(Sd_slope_df[2:4])
Sd_slope_df$resolution <- factor(Sd_slope_df$resolution, levels = c("100m resolution", "200m resolution", "300m resolution", "500m resolution", "1000m resolution"))
sd_slope_plot_trans <- ggplot(data = Sd_slope_df, 
       aes(x = mean, y = resolution, xmin = lower, xmax = upper)) +
  geom_pointrange() +
  labs(title = expression(paste("Model estimate of ", sigma, "(slope) on the estimate Shannon index")),
       x = "Coefficient Estimate",
       y = "Resolution",
       caption = "Models fit with INLA. Error bars show the 95% confidence interval.")+
  geom_vline(xintercept = as.numeric(1), col="red", alpha=0.5, lty="dashed")+
  theme_bw()

sd_slope_plot <- ggplot(data = Sd_slope_df, 
                              aes(x = mean, y = resolution, xmin = lower, xmax = upper)) +
  geom_pointrange() +
  labs(
    # title = expression(paste("Model estimate of ", sigma, "(slope) on the estimate Shannon index")),
       x = "Coefficient Estimate",
       y = "Resolution",
       caption = "Models fit with INLA using a log-link function. Error bars show the 95% confidence interval.")+
  geom_vline(xintercept = as.numeric(0), col="red", alpha=0.5, lty="dashed")+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major.y = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
sd_slope_plot


plot_grid(sd_topo_plot_trans, sd_slope_plot_trans)
plot_grid(sd_topo_plot_trans, mean_slope_plot_trans)


# glm plotting 
plot_grid(res_100_plot_pred, res_200_plot_pred, res_300_plot_pred, res_500_plot_pred, res_1000_plot_pred)


# vizualisation at 100m of the shannon and sd slop maps
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 10, "PC1278", sep="_"))
Slope <- rast("~/data/ArcDEM/sd_slope_masked_10_res.tif")

viridis_colors <- viridis::plasma(20)
na_col <- grey(0.8, alpha=0.5)
m <- ggplot() +
  geom_spatraster(data = Shannon, na.rm = TRUE, aes(fill=Shannon_10_PC1278))+ #need to change the fill variable
  # theme_map()+
  theme_minimal()+
  scale_fill_gradientn(colours = viridis_colors, na.value = na_col) +
  labs(fill = substitute(paste("Shannon index", italic("H'"))))
m  
map_div <- m+
  ggspatial::annotation_north_arrow(
    location = "bl", which_north = "true",
    pad_x = unit(0.5, "in"), pad_y = unit(0.65, "in"),
    height = unit(0.7, "cm"),
    width = unit(0.7, "cm")
  )+
  ggspatial::annotation_scale(
    location = "bl", pad_x = unit(0.5, "in"),
    pad_y = unit(0.4, "in"), 
    style="ticks"
  )
map_div

viridis_colors <- viridis::plasma(20)
na_col <- grey(0.8, alpha=0.5)
sl <- ggplot() +
  geom_spatraster(data = Slope, na.rm = TRUE, aes(fill=slope))+ #need to change the fill variable
  # theme_map()+
  theme_minimal()+
  # scale_fill_gradientn(colours = terrain.colors(10))+
  scale_fill_gradientn(colours = viridis_colors[15:1], na.value = na_col) +
  labs(fill = expression(paste(sigma,"(slope)")))
sl  
sl_div <- sl+
  ggspatial::annotation_north_arrow(
    location = "bl", which_north = "true",
    pad_x = unit(0.5, "in"), pad_y = unit(0.65, "in"),
    height = unit(0.7, "cm"),
    width = unit(0.7, "cm")
  )+
  ggspatial::annotation_scale(
    location = "bl", pad_x = unit(0.5, "in"),
    pad_y = unit(0.4, "in"), 
    style="ticks"
  )
sl_div

sent_view <- rast("~/data/output/sent_crop_view.tif")
plotRGB(sent_view, r=3, g=2, b=1, scale=10000, stretch="lin", smooth=T)
sent_view_scale <-sent_view/7

s <- ggplot() +
  geom_spatraster_rgb(data = sent_view_scale, interpolate=T, r=3, g=2, b=1)+
  # theme_map()+
  theme_minimal()
s
map_truecol <- s+
  ggspatial::annotation_scale(location="bl", pad_x=unit(0.5, "in"),
                              pad_y = unit(0.4, "in"), style="ticks", line_col="white", text_col="white")+
  
  ggspatial::annotation_north_arrow(location="bl", which_north=T, pad_x=unit(0, "in"),
                                    pad_y = unit(0.4, "in"), height = unit(0.6, "cm"), width=unit(0.6, "cm"))

map_truecol

cowplot::plot_grid(map_truecol, map_div, sl_div, ncol = 3, axis="r")


##########
# Interaction between sd(elevation) - sd(slope)??? 
##########
library(usdm)
Elev <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 10, "_res.tif", sep=""))
Slope_mean <- rast("~/data/ArcDEM/mean_slope_masked_10_res.tif")
Slope_sd <- rast("~/data/ArcDEM/sd_slope_masked_10_res.tif")
Elev_slope_mean <- c(Elev, Slope_mean)
Elev_slope_sd <- c(Elev, Slope_sd)
vif(Elev_slope_mean)
vif(Elev_slope_sd)

Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 10, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Elev <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 10, "_res.tif", sep=""))
Elev_matrix <- as.matrix(Elev, wide=T)
Slope_sd <- rast("~/data/ArcDEM/sd_slope_masked_10_res.tif")
Slope_matrix <- as.matrix(Slope_sd, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
sl <- inla.matrix2vector(Slope_matrix)
e <- inla.matrix2vector(Elev_matrix)
table(is.na(Slope_matrix))
table(is.na(Shannon_matrix))
table(is.na(Elev_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, topo=e, slope_sd = sl, node=node)
data <- data[-c(which(is.na(s))),]
hist(data$slope_sd)
hist(data$topo)
formula= Shannon_index ~ 1 + slope_sd + topo + slope_sd:topo +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_ele_sd_slope_res_100 <- inla(formula,
                                        family = "gamma",
                                        data = data,
                                        control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_ele_sd_slope_res_100)
# the WAIC and DIC is a bit higher than the simple model => maybe better to just keep the model separate.  


library(spatialEco)
corr_topo_sd_slope <- raster.modified.ttest(Elev, Slope_sd)
viridis_colors <- viridis::plasma(20)
plot(corr_topo_sd_slope$corr, type="interval", breaks=c(-1,-0.5,-0.1,0.1,0.5, 1), main="correlation between sd(elevation) and sd(slope)", col=c(viridis_colors[c(4,8,12,16)], "#009200"), plg=list(title="correlation coefficient"))



############
# Model with 10m resolution elevation data
############

# Slope
# 100m
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 10, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast("~/data/ArcDEM/sd_slope_agg10_masked_10_res.tif")
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
Gamma_shannon_sd_slope_res_100 <- inla(formula,     
                                       family = "gamma",
                                       data = data,
                                       control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_sd_slope_res_100)
save(Gamma_shannon_sd_slope_res_100, file="~/data/output/INLA_modelling/Gamma_shannon_sd_slope_agg10_res_100.Rdata")
get(load("~/data/output/INLA_modelling/Gamma_shannon_sd_slope_agg10_res_100.Rdata"))


# 200m
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 20, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast("~/data/ArcDEM/sd_slope_agg10_masked_20_res.tif")
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
Gamma_shannon_sd_slope_res_200 <- inla(formula,     
                                       family = "gamma",
                                       data = data,
                                       control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_sd_slope_res_200)
save(Gamma_shannon_sd_slope_res_200, file="~/data/output/INLA_modelling/Gamma_shannon_sd_slope_agg10_res_200.Rdata")
get(load("~/data/output/INLA_modelling/Gamma_shannon_sd_slope_agg10_res_200.Rdata"))


# 300m
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 30, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast("~/data/ArcDEM/sd_slope_agg10_masked_30_res.tif")
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
Gamma_shannon_sd_slope_res_300 <- inla(formula,     
                                       family = "gamma",
                                       data = data,
                                       control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_sd_slope_res_300)
save(Gamma_shannon_sd_slope_res_300, file="~/data/output/INLA_modelling/Gamma_shannon_sd_slope_agg10_res_300.Rdata")
get(load("~/data/output/INLA_modelling/Gamma_shannon_sd_slope_agg10_res_300.Rdata"))

# 500m
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 50, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast("~/data/ArcDEM/sd_slope_agg10_masked_50_res.tif")
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
Gamma_shannon_sd_slope_res_500 <- inla(formula,     
                                       family = "gamma",
                                       data = data,
                                       control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_sd_slope_res_500)
save(Gamma_shannon_sd_slope_res_500, file="~/data/output/INLA_modelling/Gamma_shannon_sd_slope_agg10_res_500.Rdata")
get(load("~/data/output/INLA_modelling/Gamma_shannon_sd_slope_agg10_res_500.Rdata"))


# 1000m
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 100, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
Slope <- rast("~/data/ArcDEM/sd_slope_agg10_masked_100_res.tif")
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
Gamma_shannon_sd_slope_res_1000 <- inla(formula,     
                                       family = "gamma",
                                       data = data,
                                       control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_sd_slope_res_1000)
save(Gamma_shannon_sd_slope_res_1000, file="~/data/output/INLA_modelling/Gamma_shannon_sd_slope_agg10_res_1000.Rdata")
get(load("~/data/output/INLA_modelling/Gamma_shannon_sd_slope_agg10_res_1000.Rdata"))


############
# Model with distance to water
############
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 10, "PC1278", sep="_"))
Shannon_matrix <- as.matrix(Shannon, wide=T)
dist <- rast("~/data/biodivmapR_sent/dist_to_water.tif")
dist_matrix <- as.matrix(dist, wide=T)
nrow= nrow(Shannon_matrix)
ncol= ncol(Shannon_matrix)
n = nrow*ncol
s = inla.matrix2vector(Shannon_matrix)
table(is.na(s))
s[!is.na(s)] <- s[!is.na(s)]+ 1
d <- inla.matrix2vector(dist_matrix)
table(is.na(dist_matrix))
table(is.na(Shannon_matrix))
node = 1:n
data <- data.frame(Shannon_index = s, dist = d, node=node)
data <- data[-c(which(is.na(s))),]
hist(data$dist)
plot(data$Shannon_index~data$dist)
formula= Shannon_index ~ 1 + dist +
  f(node, model="matern2d", nrow=nrow, ncol=ncol)
Gamma_shannon_dist_res_100 <- inla(formula,     
                                    family = "gamma",
                                    data = data,
                                    control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)
summary(Gamma_shannon_dist_res_100)
save(Gamma_shannon_dist_res_100, file="~/data/output/INLA_modelling/Gamma_shannon_dist_res_100.Rdata")
get(load("~/data/output/INLA_modelling/Gamma_shannon_dist_res_100.Rdata"))














############
# Plot in the thesis results 
###########
get(load("~/scratch/INLA_modelling/Gamma_shannon_sd_slope_res_100.Rdata"))
get(load("~/scratch/INLA_modelling/Gamma_shannon_sd_slope_res_200.Rdata"))
get(load("~/scratch/INLA_modelling/Gamma_shannon_sd_slope_res_300.Rdata"))
get(load("~/scratch/INLA_modelling/Gamma_shannon_sd_slope_res_500.Rdata"))
get(load("~/scratch/INLA_modelling/Gamma_shannon_sd_slope_res_1000.Rdata"))

Sd_slope_df <- rbind(Gamma_shannon_sd_slope_res_100$summary.fixed[2,c(1,3,5)],
                     Gamma_shannon_sd_slope_res_200$summary.fixed[2,c(1,3,5)],
                     Gamma_shannon_sd_slope_res_300$summary.fixed[2,c(1,3,5)],
                     Gamma_shannon_sd_slope_res_500$summary.fixed[2,c(1,3,5)],
                     Gamma_shannon_sd_slope_res_1000$summary.fixed[2,c(1,3,5)])
est <- c("100 m", "200 m", "300 m", "500 m", "1000 m")
Sd_slope_df <- cbind(as.factor(est), Sd_slope_df)
colnames(Sd_slope_df)[c(1,3:4)] <- c("resolution","lower", "upper")
# Sd_slope_df[2:4] <- exp(Sd_slope_df[2:4])
Sd_slope_df$resolution <- factor(Sd_slope_df$resolution, levels = c("100 m", "200 m", "300 m", "500 m", "1000 m"))
sd_slope_plot <- ggplot(data = Sd_slope_df, 
                        aes(x = mean, y = resolution, xmin = lower, xmax = upper)) +
  geom_pointrange(fatten=2.5) +
  labs(
    # title = expression(paste("Model estimate of ", sigma, "(slope) on the estimate Shannon index")),
    x = expression(paste("Coefficient estimates of ", sigma, " (slope)")),
    y = ""
    # ,caption = "Models fit with INLA using a log-link function. Error bars show the 95% confidence interval."
    )+
  geom_vline(xintercept = as.numeric(0), col="red", alpha=0.5, lty="dashed")+
  theme_bw() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.y = element_text(size=12), axis.text.x = element_text(size=12)
        ,panel.grid.major.x = element_blank(),  panel.grid.major.y = element_blank(), axis.title.x = element_text(vjust=-2, size=14)
        )
sd_slope_plot

# LOAD THE GOOD DATA
ggplot <- data.frame(observed= data$Shannon_index, fitted = Gamma_shannon_sd_slope_res_100$summary.fitted.values$mean[1:length(observed)])
ggplot_sd_slope_res100 <-  ggplot()+
  geom_point(data=ggplot, aes(x=fitted, y=observed), shape=21)+
  geom_abline(intercept = 0, slope = 1)+
  scale_x_continuous(breaks = c(1,2,3))+
  theme_cowplot()
ggplot_sd_slope_res100

##### coeff var 
get(load("~/scratch/INLA_modelling/Gamma_shannon_coeffvar_res_100.Rdata"))
get(load("~/scratch/INLA_modelling/Gamma_shannon_coeffvar_res_200.Rdata"))
get(load("~/scratch/INLA_modelling/Gamma_shannon_coeffvar_res_300.Rdata"))
get(load("~/scratch/INLA_modelling/Gamma_shannon_coeffvar_res_500.Rdata"))
get(load("~/scratch/INLA_modelling/Gamma_shannon_coeffvar_res_1000.Rdata"))

CV_df <- rbind(Gamma_shannon_coeffvar_res_100$summary.fixed[2,c(1,3,5)],
               Gamma_shannon_coeffvar_res_200$summary.fixed[2,c(1,3,5)],
               Gamma_shannon_coeffvar_res_300$summary.fixed[2,c(1,3,5)],
               Gamma_shannon_coeffvar_res_500$summary.fixed[2,c(1,3,5)],
               Gamma_shannon_coeffvar_res_1000$summary.fixed[2,c(1,3,5)])
est <- c("100 m", "200 m", "300 m", "500 m", "1000 m")
CV_df <- cbind(as.factor(est), CV_df)
colnames(CV_df)[c(1,3:4)] <- c("resolution","lower", "upper")
# Slope_df[,2:4] <- exp(Slope_df[,2:4])
CV_df$resolution <- factor(CV_df$resolution, levels = c("100 m", "200 m", "300 m", "500 m", "1000 m"))

CV_coeff_plot <- ggplot(data = CV_df, 
                        aes(x = mean, y = resolution, xmin = lower, xmax = upper)) +
  geom_pointrange(fatten=2.5) +
  labs(
    # title = "Model Estimates of the coefficient of variation in elevation on Shannon index estimates",
    x = expression(paste("Coefficient estimates of ", italic("CV"),"(elevation)")),
    y = ""
    # ,caption = "Models fit with INLA using a log-link function. Error bars show the 95% confidence interval."
  )+
  geom_vline(xintercept = as.numeric(0), linetype='longdash', col="red", alpha=0.5)+
  coord_cartesian(xlim = c(-0.01, 0.045))+
  theme_bw() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),  axis.text.y = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(size=12) 
        ,panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), axis.title.x = element_text(vjust=-2, size=14)
        )
CV_coeff_plot

library(grid)
library(gridExtra)
plotgrid_sd_slope_CV_ele <- plot_grid(sd_slope_plot, CV_coeff_plot, ncol=2, rel_widths = c(1.1,1))
x.grob <- textGrob("")
y.grob <- textGrob("Resolution", rot=90, gp=gpar(fontsize=14)) 
final_plot <- grid.arrange(arrangeGrob(plotgrid_sd_slope_CV_ele, left = y.grob, bottom=x.grob ))
save_plot(filename="~/data/output/final_plot/plot_grid_sd_slope_CV_ele.png", final_plot, ncol=2, base_height = 4, base_width = 5)
save_plot(filename="~/data/output/final_plot/plot_grid_sd_slope_CV_ele_with_vgrid.png", final_plot, ncol=2)


ggplot <- data.frame(observed= data$Shannon_index, fitted = Gamma_shannon_coeffvar_res_100$summary.fitted.values$mean[1:length(observed)])
ggplot_CV_res100 <-  ggplot()+
  geom_point(data=ggplot, aes(x=fitted, y=observed), shape=21)+
  geom_abline(intercept = 0, slope = 1)+
  scale_x_continuous(breaks = c(1,2,3))+
  theme_cowplot()
ggplot_CV_res100

qqplot_sd_slope_CV <- plot_grid(ggplot_sd_slope_res100, ggplot_CV_res100)
save_plot(qqplot_sd_slope_CV, filename = "~/data/output/final_plot/qqplot_sd_slope_CVele.png", base_height =4, base_width = 8, bg="white")

#### elevation
get(load("~/scratch/INLA_modelling/Gamma_shannon_topo_res100.Rdata"))
get(load("~/scratch/INLA_modelling/Gamma_shannon_topo_res200.Rdata"))
get(load("~/scratch/INLA_modelling/Gamma_shannon_topo_res300.Rdata"))
get(load("~/scratch/INLA_modelling/Gamma_shannon_topo_res500.Rdata"))
get(load("~/scratch/INLA_modelling/Gamma_shannon_topo_res1000.Rdata"))
Sd_ele_df <- rbind(Gamma_shannon_topo_res100$summary.fixed[2,c(1,3,5)],
                   Gamma_shannon_topo_res200$summary.fixed[2,c(1,3,5)],
                   Gamma_shannon_topo_res300$summary.fixed[2,c(1,3,5)],
                   Gamma_shannon_topo_res500$summary.fixed[2,c(1,3,5)],
                   Gamma_shannon_topo_res1000$summary.fixed[2,c(1,3,5)])
est <- c("100 m", "200 m", "300 m", "500 m", "1000 m")
Sd_ele_df <- cbind(as.factor(est), Sd_ele_df)
colnames(Sd_ele_df)[c(1,3:4)] <- c("resolution","lower", "upper")
Sd_ele_df$resolution <- factor(Sd_ele_df$resolution, levels = c("100 m", "200 m", "300 m", "500 m", "1000 m"))
sd_ele_plot <- ggplot(data = Sd_ele_df, 
                        aes(x = mean, y = resolution, xmin = lower, xmax = upper)) +
  geom_pointrange(fatten=2.5) +
  labs(
    # title = expression(paste("Model estimate of ", sigma, "(slope) on the estimate Shannon index")),
    x = "Coefficient estimates of mean elevation",
    y = "Resolution"
    # ,caption = "Models fit with INLA using a log-link function. Error bars show the 95% confidence interval."
  )+
  geom_vline(xintercept = as.numeric(0), col="red", alpha=0.5, lty="dashed")+
  theme_bw() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.y = element_text(size=12), axis.text.x = element_text(size=12)
        ,panel.grid.major.x = element_blank(),  panel.grid.major.y = element_blank(), axis.title.x = element_text(vjust=-2, size=14), axis.title.y = element_text(size=14),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  )
sd_ele_plot
save_plot(sd_ele_plot, filename="~/data/output/final_plot/coeff_plot_sd_ele_sent.png", base_height = 4, base_width = 5,bg="white")
save_plot(sd_ele_plot, filename="~/data/output/final_plot/coeff_plot_sd_ele_sent.svg", base_height = 4, bg="white")


Elev <- rast(paste("~/data/ArcDEM/ArcDEM_masked_", 10, "_res.tif", sep=""))
Elev_matrix <- as.matrix(Elev, wide=T)
Mean_ele <- rast(paste("~/data/ArcDEM/mean_ele_masked_", 10, "_res.tif", sep=""))
Mean_ele_matrix <- as.matrix(Mean_ele, wide=T)
e <- inla.matrix2vector(Elev_matrix)
me <- inla.matrix2vector(Mean_ele_matrix)
data <- data.frame(topo=e, mean_ele = me)
data <- data[-c(which(is.na(e))),]
data$coeff_var <- data$topo/data$mean_ele
mean_ele_sd_ele <- ggplot(data, aes(x=mean_ele, y=topo))+
  geom_point(shape=21)+
  theme_cowplot()+
  labs(x="Mean elevation", y=expression(paste(sigma, " (elevation)")))+
  theme(axis.line=element_line(linewidth =0.5), axis.text.y = element_text(size=12), axis.text.x = element_text(size=12),
        axis.title.x = element_text(vjust=-2,size=14),axis.title.y=element_text(size=14), 
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
        )
mean_ele_sd_ele
save_plot(mean_ele_sd_ele, filename = "~/data/output/final_plot/scatterplot_mean_ele_sd_ele.png", bg="white", base_height = 4, base_asp = 1.4)

plotgrid_sd_ele_scatterplot <- plot_grid(sd_ele_plot, mean_ele_sd_ele, ncol=2, rel_widths = c(1,1))
save_plot(plotgrid_sd_ele_scatterplot, filename = "~/data/output/final_plot/plotgrid_coeff_sd_ele_scatterplot.png", bg="white", ncol=2, base_height = 4, base_width = 5)

