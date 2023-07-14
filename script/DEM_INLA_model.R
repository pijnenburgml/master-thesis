library(sp)  # package to work with spatial data
library(rgdal)  
require(spdep)

# install.packages("INLA",
#                  repos = c(getOption("repos"),
#                            INLA = "https://inla.r-inla-download.org/R/stable"),
#                  dep = T)
library(INLA)
library(terra)
Shannon_map <- rast("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi/SPCA/ALPHA/Shannon_10")
table(is.na(Shannon_map[1:120540]))
Sentinel_map_poly <- as.polygons(Shannon_map, round=F, aggregate=F, extent=F)
# writeVector(Sentinel_map_poly, filename = "~/data/output/Sentinel_lattice.shp")

Sentinel_sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID = "1")))
Sentinel_sf_poly <- st_read("~/data/output/Sentinel_lattice.shp")

Sentinel_lattice <- readOGR("~/data/output/model_object.shp")
Sentinel_data <- Sentinel_lattice@data
plot(Sentinel_data$Shannon_10, Sentinel_data$sd_topo)
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


formula <- Shannon_10 ~ 1 + sd_topo + # fixed effect
  f(spatial, model = "bym",       # spatial effect: ZONE_CODE is a numeric identifier for each area in the lattice  (does not work with factors)
    graph = Lattice.adj)

hist(Sentinel_lattice$Shannon_10)
sample <- sample(Sentinel_lattice$Shannon_10, size=5000)
hist(sample)
shapiro.test(sample)
qqnorm(Sentinel_lattice$Shannon_10)
qqline(Sentinel_data$Shannon_10)

Mod_Lattice <- inla(formula,     
                    family = "gaussian", # have to change the family, not gaussian
                    data = Sentinel_data,
                    control.compute = list(cpo = T, dic = T, waic = T))  
# CPO, DIC and WAIC metric values can all be computed by specifying that in the control.compute option
# These values can then be used for model selection purposes if you wanted to do that

# Check out the model summary
summary(Mod_Lattice)
