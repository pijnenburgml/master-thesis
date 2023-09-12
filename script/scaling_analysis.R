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

Gaussian_loglink_shannon_topo <- inla(formula,family = "gaussian",
                          control.family=list(link='identity'),
                          data = data,
                          control.compute = list(cpo = T, dic = T, waic = T, return.marginals.predictor=TRUE), verbose=TRUE)

summary(Gaussian_loglink_shannon_topo)
# save(Gaussian_loglink, file="~/data/output/INLA_modelling/Gaussian_model_loglink_shannon_log_topo.Rdata")
Gaussian_loglink <- get(load("~/data/output/INLA_modelling/Gaussian_model_loglink_shannon_log_topo.Rdata"))


plot_inla_residuals(Gaussian_loglink_withNA, observed=data$Shannon_index)

ggplot_inla_residuals(Model_Lattice_no_ndvi, observed=observed)
ggplot_inla_residuals2(Model_Lattice_no_ndvi, observed, se = FALSE)


formula= Shannon_index ~ 1+ log(topo)+
  f(node, model="matern2d", nrow=nrow, ncol=ncol)

Gaussian_shannon_log_topo <- inla(formula,family = "gaussian",
                                      control.family=list(link='identity'),
                                      data = data,
                                      control.compute = list(cpo = T, dic = T, waic = T, return.marginals.predictor=TRUE), verbose=TRUE)


summary(Gaussian_shannon_log_topo)
plot(data$Shannon_index~log(data$topo))
#####
# Gamma
#####

# seems better than gaussian

formula= Shannon_index ~ 1+ log(topo)+
  f(node, model="matern2d", nrow=nrow, ncol=ncol)

Gamma_shannon_log_topo <- inla(formula,     
                            family = "gamma", # have to change the family, not gaussian, try gamma
                            data = data,
                            control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)

summary(Gamma_shannon_log_topo)
# save(Gamma_shannon_log_topo, file="~/data/output/INLA_modelling/Gamma_shannon_log_topo.Rdata")
Gamma_shannon_log_topo <- get(load("~/data/output/INLA_modelling/Gaussian_model_loglink_shannon_log_topo.Rdata"))

formula= Shannon_index ~ 1+ topo+
  f(node, model="matern2d", nrow=nrow, ncol=ncol)

Gamma_shannon_topo <- inla(formula,     
                               family = "gamma",
                               data = data,
                               control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE)

summary(Gamma_shannon_topo)

exp(Gamma_shannon_topo$summary.fixed)
Gamma_shannon_topo_res100 <- Gamma_shannon_topo

save(Gamma_shannon_topo_res100, file="~/data/output/INLA_modelling/Gamma_shannon_topo_res100.Rdata")
get(load("~/data/output/INLA_modelling/Gamma_shannon_topo_res100.Rdata"))

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

save(Gamma_shannon_topo_res200, file="~/data/output/INLA_modelling/Gamma_shannon_topo_res200.Rdata")
get(load("~/data/output/INLA_modelling/Gamma_shannon_topo_res200.Rdata"))

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

save(Gamma_shannon_topo_res300, file="~/data/output/INLA_modelling/Gamma_shannon_topo_res300.Rdata")
get(load("~/data/output/INLA_modelling/Gamma_shannon_topo_res300.Rdata"))

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
get(load("~/data/output/INLA_modelling/Gamma_shannon_topo_res500.Rdata"))


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
get(load("~/data/output/INLA_modelling/Gamma_shannon_topo_res1000.Rdata"))

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
  xlab("coefficient sd(topography)") +
  ylab("Estimate")+
  coord_flip()+
  theme_bw()+
  geom_vline(xintercept = 0, linetype='longdash', col="red")
sd_topo_plot_trans

par(mfrow=c(1,2))
sd_topo_plot
sd_topo_plot_trans
plot_grid(sd_topo_plot, sd_topo_plot_trans)



##################
# simple glm regression
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


