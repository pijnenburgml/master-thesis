library(sp)  # package to work with spatial data
library(rgdal)  
require(spdep)

# install.packages("INLA",
#                  repos = c(getOption("repos"),
#                            INLA = "https://inla.r-inla-download.org/R/stable"),
#                  dep = T)
library(INLA)
library(terra)
library(ggplot2)
library(ggregplot)
library(tidyverse)
library(RColorBrewer)
library(INLAutils)


########
# NDVI for modelilng
########
SentMask <- rast("~/data/biodivmapR_sent/mask_sent2_final_NA")
NDVI <- rast("~/data/output/Sentinel_NDVI.tif")
tile_DEM_crop <- rast("~/data/ArcDEM/tile_DEM_crop.tif")
tile_DEM_100m <- aggregate(tile_DEM_crop, fact=100/2, fun="sd") #1.993736
Shannon_map <- rast("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi/SPCA/ALPHA/Shannon_10")
Shannon_proj_DEM <- terra::project(Shannon_map, tile_DEM_100m, method="bilinear")
Shannon_masked <- mask(Shannon_proj_DEM, tile_DEM_100m_masked)
tile_DEM_mask_resample <- terra::project(tile_DEM_crop, SentMask, method="bilinear")

NDVI_100b100 <- resample(NDVI, Shannon_masked, method="bilinear")
plot(NDVI_100b100)
NDVI_100b100_masked <- mask(NDVI_100b100, Shannon_masked)
plot(NDVI_100b100_masked)
plot(Shannon_masked, add=T, col="grey")
plot(tile_DEM_100m_masked, add=T, col="red")
Sentinel_lattice <-readOGR("~/data/ArcDEM/model_object.shp")

############
# INLA modelling
############
setwd("~/data/output/INLA_modelling/")
Shannon_map <- rast("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_10")
table(is.na(Shannon_map[1:120540]))
plot(Shannon_map)
# Sentinel_map_poly <- as.polygons(Shannon_map, round=F, aggregate=F, extent=F)
# writeVector(Sentinel_map_poly, filename = "~/data/output/Sentinel_lattice.shp")

# Sentinel_sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID = "1")))
# Sentinel_sf_poly <- st_read("~/data/output/Sentinel_lattice.shp")

Sentinel_lattice <- readOGR("~/data/ArcDEM/model_object.shp")
Sentinel_lattice$NDVI <- unlist(as.data.frame(NDVI_100b100_masked)[1])
# Sentinel_lattice <- readOGR("~/data/output/model_object_no_0.shp")
Sentinel_data <- Sentinel_lattice@data
plot(x=Sentinel_data$Shannon_10, y=Sentinel_data$NDVI)
plot(Sentinel_data$Shannon_10, Sentinel_data$sd_topo)
plot(Sentinel_data$Shannon_10, log(Sentinel_data$sd_topo))

lattice_temp <- poly2nb(Sentinel_lattice)
# listw <- nb2listw(lattice_temp, zero.policy = T)
# mat <- nb2mat(lattice_temp, zero.policy = T)
nb2INLA("Lattice.graph", lattice_temp) # create the adjacency matrix in INLA format
Lattice.adj <- paste(getwd(),"/Lattice.graph",sep="") # name the object
inla.setOption(scale.model.default = F)
H <- inla.read.graph(filename = "Lattice.graph")  # and save it as a graph
image(inla.graph2matrix(H), xlab = "", ylab = "")
# don't really get what's happening....

mat_inla <- inla.graph2matrix(H)
image(inla.graph2matrix(H)[1:500, 1:500])

# Sentinel_data_point_matrix <- extract(x=Shannon_map, y=Sentinel_map_poly)
# spplot(obj = Sentinel_lattice, zcol = "Shannon_10")

formula <- Shannon_10 ~ 1 + log(sd_topo) + # fixed effect
  f(spatial, model = "bym",       # spatial effect: ZONE_CODE is a numeric identifier for each area in the lattice  (does not work with factors)
    graph = Lattice.adj)

hist(Sentinel_lattice$Shannon_10, breaks=seq(0, 3, by=0.1))
hist(Sentinel_data$sd_topo, breaks=seq(0,19, by=0.1))
hist(log(Sentinel_data$sd_topo), breaks=seq(-4, 3, by=0.1))
hist(Sentinel_data$NDVI) # What to do with this distribution, not Gaussian => not good to put in a model with gaussian distribution 
Shannon_10_sample <- sample(Sentinel_lattice$Shannon_10, size=5000)
hist(Shannon_10_sample)
shapiro.test(Shannon_10_sample)
qqnorm(Sentinel_data$Shannon_10)
qqline(Sentinel_data$Shannon_10)

sd_topo_sample <- sample(Sentinel_data$sd_topo, size=5000)
hist(sd_topo_sample)
mean(sd_topo_sample)
var(sd_topo_sample)
plot(Shannon_10_sample~log(sd_topo_sample))

set.seed(10)
random_sample <- sample(1:nrow(Sentinel_data), size=5000)
Sentinel_data_sample <- Sentinel_data[random_sample,]
hist(Sentinel_data_sample$Shannon_10)
hist(Sentinel_data_sample$sd_topo)
mean(Sentinel_data_sample$Shannon_10)
mean(Sentinel_data_sample$sd_topo)
var(Sentinel_data_sample$Shannon_10)
var(Sentinel_data_sample$sd_topo)

mean(Sentinel_data$Shannon_10)
mean(Sentinel_data$sd_topo)
var(Sentinel_data$Shannon_10)
var(Sentinel_data$sd_topo)
plot(Sentinel_data_sample$Shannon_10~log(Sentinel_data_sample$sd_topo))
Sentinel_data_sample$log_sd_topo <- log(Sentinel_data_sample$sd_topo)
ggplot(Sentinel_data_sample, aes(log_sd_topo,Shannon_10)) + 
  # Add points to the plot
  geom_point() +  
  geom_smooth(se=T,
              linewidth=1, color="red",
              linetype = "dashed")



# colnames(Sentinel_data)[1] <- "y"
formula <- Shannon_10 ~ 1 + log(sd_topo) + NDVI +# fixed effect
  f(spatial, model = "bym",       # spatial effect: ZONE_CODE is a numeric identifier for each area in the lattice  (does not work with factors)
    graph = Lattice.adj)

Model_Lattice <- inla(formula,     
                    family = "gaussian", # have to change the family, not gaussian, try gamma
                    data = Sentinel_data,
                    control.compute = list(cpo = T, dic = T, waic = T, return.marginals.predictor=TRUE), verbose=TRUE, scale=1)

formula <- Shannon_10 ~ 1 + log(sd_topo)+ # fixed effect
  f(spatial, model = "bym",       # spatial effect: ZONE_CODE is a numeric identifier for each area in the lattice  (does not work with factors)
    graph = Lattice.adj)

Model_Lattice_no_ndvi <- inla(formula,     
                    family = "gaussian", # have to change the family, not gaussian, try gamma
                    data = Sentinel_data,
                    control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE, scale=1)

summary(Model_Lattice_no_ndvi)
# problem with 0 in the dataset with gamma ditribution family
save(Model_Lattice_no_ndvi, file="~/data/output/INLA_modelling/gaussian_model_shannon_log_topo.Rdata")
Model_Lattice_no_ndvi <- get(load("~/data/output/INLA_modelling/gaussian_model_shannon_log_topo.Rdata"))

nrow(Model_Lattice_no_ndvi$summary.fixed)
coef_df <- data.frame()
for(i in 1:nrow(Model_Lattice_no_ndvi$summary.fixed)){
  name <- NULL
  fixed <- NULL
  d <- NULL
  name <- Model_Lattice_no_ndvi$summary.fixed[i,] %>% rownames()
  fixed <- Model_Lattice_no_ndvi$summary.fixed[i,] %>% as.numeric()
  var_df <- data.frame(
    var = name,
    lower = fixed[3],
    upper = fixed[5],
    estimate = fixed[1]
  )
  
  coef_df <- rbind(coef_df, var_df)
}
name <- "sd_topo"
fixed <- exp(Model_Lattice_no_ndvi$summary.fixed[2,]) %>% as.numeric()
var_df <- data.frame(
  var = name,
  lower = fixed[3],
  upper = fixed[5],
  estimate = fixed[1]
)
coef_df <- rbind(coef_df, var_df)

ggplot(data=coef_df, aes(estimate, var)) +
    geom_point()+
    # geom_linerange(aes(xmin = lower, xmax = upper))+
    geom_vline(xintercept = 0, lty = 2, col="red")

ggplot(data=coef_df, aes(estimate, var))+
  geom_pointrange(aes(xmin = lower, xmax = upper))+
  geom_vline(xintercept = 0, lty = 2, col="red")

# CPO, DIC and WAIC metric values can all be computed by specifying that in the control.compute option
# These values can then be used for model selection purposes if you wanted to do that

plot(Model_Lattice_no_ndvi, plot.fixed.effects=FALSE, plot.lincomb=FALSE, plot.random.effects=FALSE,
     plot.hyperparameters=FALSE, plot.predictor=FALSE, plot.q=FALSE, plot.cpo=TRUE,
     single=FALSE)

inla.show.hyperspec(Model_Lattice_no_ndvi) 

plot(Sentinel_data$Shannon_10~Sentinel_data$sd_topo)
points(Sentinel_data[which(Model_Lattice_no_ndvi$cpo$failure!=0),1]~Sentinel_data[which(Model_Lattice_no_ndvi$cpo$failure!=0),2], col="red")


plot(Model_Lattice_no_ndvi, plot.fixed.effects=FALSE, plot.lincomb=FALSE, plot.random.effects=FALSE,
     plot.hyperparameters=FALSE, plot.predictor=TRUE, plot.q=FALSE, plot.cpo=FALSE,
     single=FALSE)


# Check out the model summary
plot(Model_Lattice_no_ndvi)
summary(Model_Lattice_no_ndvi)
plot(Sentinel_data$Shannon_10, Model_Lattice_no_ndvi$summary.fitted.values$mean, main="Fitting result")
observed <- Sentinel_data$Shannon_10
plot_inla_residuals(Model_Lattice_no_ndvi, observed=observed)
ggplot_inla_residuals(Model_Lattice_no_ndvi, observed=observed)
ggplot_inla_residuals2(Model_Lattice_no_ndvi, observed, se = FALSE)


library(ggplot2)
intercept <- Model_Lattice_no_ndvi$marginals.fixed[[1]]
ggplot(data.frame(inla.smarginal(intercept)), aes(x, y)) +
  geom_line() +
  theme_bw()

quant <- inla.qmarginal(0.05, intercept)
quant
inla.pmarginal(quant, intercept)
ggplot(data.frame(inla.smarginal(intercept)), aes(x, y)) +
  geom_line() +
  geom_area(data = subset(data.frame(inla.smarginal(intercept)),
                          x < quant),
            fill = "black") +
  theme_bw()



topo <- exp(Model_Lattice_no_ndvi$marginals.fixed[[2]])
ggplot(data.frame(inla.smarginal(topo)), aes(x, y)) +
  geom_line() +
  theme_bw()

marg.topo.trans <- inla.tmarginal(function(x) exp(x),
                                Model_Lattice_no_ndvi$marginals.fixed[[2]])
ggplot(data.frame(inla.smarginal(marg.topo.trans)), aes(x, y)) +
  geom_line() +
  theme_bw()

library(remotes)
remotes::install_github("gfalbery/ggregplot")

library(ggregplot)
library(dplyr)
Efxplot(Model_Lattice_no_ndvi)
# error doesn't find the function %<>%



################
# model with Gamma family and log-link function 
################

Sentinel_lattice <- readOGR("~/data/ArcDEM/model_object.shp")
Sentinel_lattice$Shannon_10 <- Sentinel_lattice$Shannon_10+1
Sentinel_data <- Sentinel_lattice@data
which(Sentinel_data$Shannon_10 == 0)
min(Sentinel_data$Shannon_10)
hist(Sentinel_data$Shannon_10)
plot(Sentinel_data$Shannon_10, Sentinel_data$sd_topo)
plot(Sentinel_data$Shannon_10, log(Sentinel_data$sd_topo))
lattice_temp <- poly2nb(Sentinel_lattice)
# listw <- nb2listw(lattice_temp, zero.policy = T)
# mat <- nb2mat(lattice_temp, zero.policy = T)
nb2INLA("Lattice.graph", lattice_temp) # create the adjacency matrix in INLA format
Lattice.adj <- paste(getwd(),"/Lattice.graph",sep="") # name the object
inla.setOption(scale.model.default = F)
H <- inla.read.graph(filename = "Lattice.graph")  # and save it as a graph
image(inla.graph2matrix(H), xlab = "", ylab = "")
formula <- Shannon_10 ~ 1 + log(sd_topo) + # fixed effect
  f(spatial, model = "bym",       # spatial effect: ZONE_CODE is a numeric identifier for each area in the lattice  (does not work with factors)
    graph = Lattice.adj)

model_lattice_gamma <- inla(formula,     
                              family = "gamma", # have to change the family, not gaussian, try gamma
                              data = Sentinel_data,
                              control.compute = list(cpo = T, dic = T, waic = T,return.marginals.predictor=TRUE), verbose=TRUE, scale=1)

save(model_lattice_gamma, file="~/data/output/INLA_modelling/gamma_model_shannon_log_topo.Rdata")
get(load("~/data/output/INLA_modelling/gamma_model_shannon_log_topo.Rdata"))

summary(model_lattice_gamma)
observed <- Sentinel_data$Shannon_10
plot_inla_residuals(model_lattice_gamma, observed=observed)


################
# test for model choice
################

reg <- lm(Shannon_10~log(sd_topo), data=Sentinel_data[,1:2])
par(mfrow=c(1,2)) 
plot(reg, which=c(1,3), pch=20) 
summary(lm(log(resid(reg)^2)~log(fitted(reg)))) 
# log(fitted(reg))  -52.449      4.459  -11.76   <2e-16 ***
# wtf is that relationship


#################
# model with gaussian family and log-link function 
#################
Sentinel_lattice <- readOGR("~/data/ArcDEM/model_object.shp")
Sentinel_lattice$Shannon_10 <- Sentinel_lattice$Shannon_10+1
Sentinel_data <- Sentinel_lattice@data
lattice_temp <- poly2nb(Sentinel_lattice)
nb2INLA("Lattice.graph", lattice_temp) # create the adjacency matrix in INLA format
Lattice.adj <- paste(getwd(),"/Lattice.graph",sep="") # name the object
inla.setOption(scale.model.default = F)
H <- inla.read.graph(filename = "Lattice.graph")  # and save it as a graph

formula <- Shannon_10 ~ 1 + log(sd_topo) + # fixed effect
  f(spatial, model = "bym",       # spatial effect: ZONE_CODE is a numeric identifier for each area in the lattice  (does not work with factors)
    graph = Lattice.adj)


model_lattice_gaussian_loglink <- inla(formula,     
                      family = "gaussian",
                      control.family=list(link='log'),
                      data = Sentinel_data,
                      control.compute = list(cpo = T, dic = T, waic = T, return.marginals.predictor=TRUE), verbose=TRUE, scale=1)

summary(model_lattice_gaussian_loglink)
# save(model_lattice_gaussian_loglink, file="~/data/output/INLA_modelling/gaussian_model_loglink_shannon_log_topo.Rdata")
Model_lattice_gaussian_loglink <- get(load("~/data/output/INLA_modelling/gaussian_model_loglink_shannon_log_topo.Rdata"))

summary(model_lattice_gaussian_loglink)
observed <- Sentinel_data$Shannon_10
plot_inla_residuals(model_lattice_gaussian_loglink, observed=observed)














