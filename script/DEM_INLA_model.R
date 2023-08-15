library(sp)  # package to work with spatial data
library(rgdal)  
require(spdep)

# install.packages("INLA",
#                  repos = c(getOption("repos"),
#                            INLA = "https://inla.r-inla-download.org/R/stable"),
#                  dep = T)
library(INLA)
library(terra)



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

# problem with 0 in the dataset




# CPO, DIC and WAIC metric values can all be computed by specifying that in the control.compute option
# These values can then be used for model selection purposes if you wanted to do that

# Check out the model summary
summary(Model_Lattice)
plot(Sentinel_data$Shannon_10, Model_Lattice$summary.fitted.values$mean, main="Fitting result")
observed <- Sentinel_data$Shannon_10
library(INLAutils)
plot_inla_residuals(Model_Lattice, observed=observed)
residual_plot

