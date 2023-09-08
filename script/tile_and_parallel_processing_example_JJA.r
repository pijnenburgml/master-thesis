# Example script for generating and processing tile geometries
# Jakob Assmann jakob.assmann@uzh.ch 27 June 2023

# Dependencies
library(sf)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(terra)
library(tidyterra)
library(devtools)
library(gdalUtils)

setwd("~/data")
area_interest <- st_read("areas_of_interest.gpkg")
area_interest_proj <- st_transform(area_interest,32613)
sent_aoi_stack_crop <- rast("output/sent_crop_view.tif")
ext_aoi <- ext(area_interest_proj)
ext_sent <- ext(sent_aoi_stack_crop)


#########
# Target area covering the whole Aviris fligth stripe 
#########
target_area <- data.frame(x = c(473557.1, 569724.9, 569724.9, 473557.1), #this is the Sentinel-2 outline.
                          y = c(7639904, 7639904, 7700841 ,7700841)) %>%
  st_as_sf(
    coords = c("x", "y"),
    crs = 32613 # UTM Zone 12 -> so coordiantes are in metres
  ) %>% summarize() %>%
  st_convex_hull()

#########
# Target area covering the aoi --> Sentinel tile
#########
target_area <- data.frame(x = c(476480, 517470, 517470, 476480), #this is the Sentinel-2 outline.
                          y = c(7657730, 7657730,  7687110, 7687110)) %>%
    st_as_sf(
    coords = c("x", "y"),
    crs = 32613 # UTM Zone 12 -> so coordiantes are in metres
  ) %>% summarize() %>%
  st_convex_hull()

plotRGB(sent_aoi_stack_crop, r=3, g=2, b=1, scale=10000, stretch="lin")
plot(target_area, add=T, col="red")


# Visualise the polygon to check everything worked
ggplot() +
    geom_sf(data = target_area) +
    theme_map()

# Set grid size (asumming square tiles)
tile_width <- 1000 # 1 km tiles

# Arbitarily define the bottom left of the bounding box of our
# target area as the origin of the grid
grid_origin <- st_bbox(target_area)[c("xmin", "ymin")]

# Calculate number of grid cells required to cover the are of interest
n_cells_x <- ceiling((st_bbox(target_area)["xmax"] - st_bbox(target_area)["xmin"]) / tile_width)
n_cells_y <- ceiling((st_bbox(target_area)["ymax"] - st_bbox(target_area)["ymin"]) / tile_width)
# st_bbox --> return xmin, ymin, xmax, ymax => here select with[] either one of them

# wirte a quick helper function to generate a polygon for a tile
generate_tile_poly <- function(cell_id_x, cell_id_y, grid_origin, tile_width) {
    # Generate a matrix of five points that form the corners of the polygon 
    # plus the start corner to make the geometry closed
    matrix(nrow = 5, ncol = 2, byrow = TRUE,
    data = c(
            grid_origin["xmin"] + (cell_id_x - 1) * tile_width, grid_origin["ymin"] + (cell_id_y - 1) * tile_width,
            grid_origin["xmin"] + (cell_id_x) * tile_width, grid_origin["ymin"] + (cell_id_y - 1) * tile_width,
            grid_origin["xmin"] + (cell_id_x) * tile_width, grid_origin["ymin"] + (cell_id_y) * tile_width,
            grid_origin["xmin"] + (cell_id_x - 1) * tile_width, grid_origin["ymin"] + (cell_id_y) * tile_width,
            grid_origin["xmin"] + (cell_id_x - 1) * tile_width, grid_origin["ymin"] + (cell_id_y - 1) * tile_width)) %>%
            # Cast into a list
            list() %>%
                # Convert to polygon
                st_polygon() %>%
                # Convert to sf with a crs
                st_sfc(crs = 32613) %>%
                st_sf() %>%
                # add identifiers to the geometry
                mutate(
                    cell_id = paste0("x_", cell_id_x, "_y_", cell_id_y),
                    cell_id_x = cell_id_x,
                    cell_id_y = cell_id_y
                )
}

# Confirm that it works
ggplot() +
    geom_sf(data = target_area)+
    geom_sf(data = generate_tile_poly(1, 1, grid_origin, tile_width),col="red") +
    theme_map()

# Generate polygons for all tiles:

# Use expand.grid to get all combinaton of grid cell ids
grid_cell_combos <- expand.grid(x = 1:n_cells_x, y = 1:n_cells_y) %>%
    data.frame()

# Use map 2 to apply the helper function
grid_polygons <- map2(grid_cell_combos$x, grid_cell_combos$y, .f = generate_tile_poly, grid_origin = grid_origin, tile_width = tile_width) %>%
    # bind list of sfcs into one
    bind_rows()

# Plot to check it worked
# mask <- rast("~/master-thesis/sent_output/mask_sent2_final.tif")
ggplot() +
  geom_sf(
    data = grid_polygons,
    aes(fill = cell_id)
  ) +
  theme_map() +
  theme(legend.position = "none")
  # +geom_spatraster(data=mask, alpha=0.3)

#Save the tile boundaries as shapefile 
# st_write(grid_polygons, "~/data/output/tilling/grid_polygons_aoi.shp")
# st_write(grid_polygons, "~/data/output/tilling/grid_polygons_aviris.shp")
grid_polygons <- st_read("~/data/output/tilling/grid_polygons_aoi.shp")
# grid_polygons <- st_read("~/data/output/tilling/grid_polygons_aviris.shp")
##############################################################

###################
# Get the fligth stripe boundaries
###################
boundary_dir <- "~/data/Aviris_data/boundaries"
boundaries <- list.files(path=boundary_dir, pattern="KML", full.names = T)
# try <- st_read("~/data/Aviris_data/boundaries/ang20190801t143347_outline_KML.kml")
# try_proj <- st_transform(try,32613)
# st_intersects(grid_polygons[1,], try_proj$geometry)
# 
# #plot enabling to tell which tile cell overlap with a flight stripe (here 3347 = try_proj)
# ggplot(data=grid_polygons) +
#   geom_sf(
#     # aes(fill = cell_id, alpha=0.1)
#   ) +
#   geom_sf_label(aes(label=cell_id), cex=0.8)+
#   geom_sf(data=boundaries_multipolygon, aes(fill="red", alpha=0.1))
# st_intersects(grid_polygons[grep(pattern="x_1_y_17", grid_polygons$cell_id),], try_proj)

######################
# make multipolygon from all the flight strip boundaries
######################

boundaries_sf <- lapply(boundaries, FUN=st_read)
for(i in 1:length(boundaries_sf)){
  boundaries_sf[[i]] <- st_transform(boundaries_sf[[i]],32613) 
}

# Get names
# substr(boundaries, 42, 75)
boundaries_sf[[1]]$Name
boundaries_name <- c()
for(i in 1:length(boundaries_sf)){
  boundaries_name[i] <- boundaries_sf[[i]]$Name
}
names(boundaries_sf) <- boundaries_name

# get all the objects from the boundaries_sf list into a  
cat(paste(rep("boundaries_sf[[", 16), 1:16, rep("]]$geometry", 16), sep=""), sep=",")
# copy-paste output
boundaries_multipolygon <- c(boundaries_sf[[1]]$geometry, boundaries_sf[[2]]$geometry,boundaries_sf[[3]]$geometry,boundaries_sf[[4]]$geometry,
  boundaries_sf[[5]]$geometry,boundaries_sf[[6]]$geometry,boundaries_sf[[7]]$geometry,boundaries_sf[[8]]$geometry,
  boundaries_sf[[9]]$geometry,boundaries_sf[[10]]$geometry,boundaries_sf[[11]]$geometry,boundaries_sf[[12]]$geometry,
  boundaries_sf[[13]]$geometry,boundaries_sf[[14]]$geometry,boundaries_sf[[15]]$geometry,boundaries_sf[[16]]$geometry)%>% st_cast("MULTIPOLYGON")
plot(boundaries_multipolygon)
names(boundaries_multipolygon) <- boundaries_name


########
# Try to apply gdal warp to get a piece of flight strip rectified I could run the random forest on
########
bbox <- dplyr::filter(grid_polygons, cell_id=="x_35_y_17") %>% st_bbox()
source_path <-"/scratch/mpijne/reflectance_data/ang20190801t160747_rfl"
sf::gdal_utils("warp",
               source = source_path, 
               destination = paste0(source_path, "_x_35_y_17"),
               options = c(
                 # outpufile = ENVI file
                 "-of", "ENVI",
                 "-t_srs", "epsg:32613",
                 # target extent (in target crs)
                 "-te",
                 bbox[1], # xmin
                 bbox[2], # ymin
                 bbox[3], # xmax
                 bbox[4] # ymax
               )
)


###############
# Using everything for parallel processing (on li)
###############
library(pbapply)
library(parallel)
cl <- 4 # change with the number of cluster 

# then run things in parallel 
pblapply(grid_polygons$cell_id,
         function(cell) {
           file <- names(boundaries_multipolygon)[which(st_intersects(boundaries_multipolygon, dplyr::filter(grid_polygons, cell_id==cell), sparse = F)==T)]
           paths <- paste0("/scratch/mpijne/reflectance_data/", substr(file, 0, 18), "_rfl")
           bbox <- dplyr::filter(grid_polygons, cell_id==cell) %>% st_bbox()
           if (length(file)>0) { #need a loop or an apply
             for (x in paths) {
             sf::gdal_utils("warp",
                            source = x, 
                            destination = paste(x, cell, sep="_"),
                            options = c(
                              # outpufile = ENVI file
                              "-of", "ENVI",
                              "-t_srs", "epsg:32613",
                              # target extent (in target crs)
                              "-te",
                              bbox[1], # xmin
                              bbox[2], # ymin
                              bbox[3], # xmax
                              bbox[4] # ymax
                              )
                            ) 
             }
             update_paths <- paste(paths, cell, sep="_") 
             sf::gdal_utils("warp", source = update_paths, destination = paste0("/scratch/mpijne/reflectance_data/", cell),
                            options = c("-of", "ENVI", "-t_srs", "epsg:32613",  "-co", "INTERLEAVE=BIL"))
           }
          }
          # cl = 4
)

# check if trial with test file worked: 
tile1 <- rast("/scratch/mpijne/reflectance_data/ang20190801t143347_rfl_x_10_y_14")
tile2 <- rast("/scratch/mpijne/reflectance_data/ang20190801t144718_rfl_x_10_y_14")
plotRGB(tile1, r=54, g=36, b=20, stretch="lin", colNA="red")
plotRGB(tile2, r=54, g=36, b=20, stretch="lin", colNA="red")
tile_x10y14 <- rast("/scratch/mpijne/reflectance_data/x_10_y_14")
plotRGB(tile_x10y14, r=54, g=36, b=20, stretch="lin")

# comparison with mosaic not possible because different resolution
# tile_x10y14 <- mosaic(tile1, tile2)
# writeRaster(tile_x35y17, "/scratch/mpijne/reflectance_data/x_35_y_17_comparison.envi")

#######
# Masking Aviris data in R
#######

# Clouds
RF_model <- randomForest(x=training_data[,2:426], y=training_data$ID)
pix_class <- terra::predict(object=tile_x10y14, model=RF_model)

# Shade
strip_4159_crop <- rast("./Aviris_data/strip_4159_cropped_2_rect")
NIR_idx <- names(strip_4159_crop) %>% grep(pattern="^8")
NIR_average <- mean(strip_4159_crop[[NIR_idx]])
hist(NIR_average, breaks=seq(0,0.35, by=0.01))
# cut hard decision, at 0.1 can make sense but less as well, something like 0.05 could also be justifiable.  
shade_mask <- ifel(NIR_average<0.05, 0, 1)
plot(shade_mask)
# NDWI
green_idx <- names(strip_4159_crop) %>% grep(pattern="^5")
green_average <- mean(strip_4159_crop[[green_idx]])
NDWI <- (green_average-NIR_average)/(green_average+NIR_average)
plot(NDWI)
hist(NDWI)
water_mask <- ifel(NDWI>0.1, 0, 1)
plot(water_mask)
par(mfrow=c(1,2))
plot(shade_mask)
plot(water_mask)
plotRGB(strip_4159_crop, r=names(strip_4159_crop)[54], g=names(strip_4159_crop)[36], b=names(strip_4159_crop)[20], stretch="lin")

RF_model <-readRDS("~/data/output/random_forest_cloud_detection/cloud_classifier_RF.RData")
build_area_mask <- rast("~/data/output/build_are_mask.tif")
sent_aoi_stack_crop <- rast("output/sent_crop_view.tif")
temp_rast <- rast(ext(sent_aoi_stack_crop), resolution = res(tile))
build_area_mask_5res <- resample(build_area_mask, temp_rast, method="near")

pblapply(grid_polygons$cell_id,
         function(cell) {
           tile <- rast(file.path("/scratch/mpijne/reflectance_data", cell))
           NIR_average <- mean(tile[[grep(names(tile), pattern = "^8")]])
           # hist(NIR_average, breaks=seq(terra::minmax(NIR_average)[1],terra::minmax(NIR_average)[2]+0.05, by=0.01))
           shade_mask <- ifel(NIR_average<0.05, 0, 1)
           green_average <- mean(tile[[grep(names(tile), pattern="^5")]])
           NDWI <- (green_average-NIR_average)/(green_average+NIR_average)
           water_mask <- ifel(NDWI>0.1, 0, 1)
           names(tile) <- 1:425
           cloud_class <- terra::predict(object=tile, model=RF_model)
           cloud_mask <- ifel(cloud_class=="cl", 0, 1)
           build_area_mask <- rast("~/data/output/build_are_mask.tif")
           temp_rast <- rast(ext(sent_aoi_stack_crop), resolution = res(tile))
           build_area_mask_5res <- resample(build_area_mask, temp_rast, method="near")
           building_mask <- crop(build_area_mask_5res, tile)
           mask <- mosaic(shade_mask, water_mask, cloud_mask, building_mask, fun="min")
           # mask <- ifel(mask==0, NA, 1) #is it necessary?
           writeRaster(mask, filename = paste("/scratch/mpijne/reflectance_data/", cell, "_mask", sep=""),filetype="ENVI", gdal="INTERLEAVE=BSQ", overwrite=T, datatype="INT1U")
                         }
         cl = cl # the name of your cluster
)


############################
# parallel processing to create Output_30_PCs
############################
pblapply(grid_polygons$cell_id,
         function(cell) {
            
           
         }
)






###############
# Code from Jakob
################

# Quick example on how to iterate about tiles in parallel
# For this I prefer using the library plapply
library(pbapply)
# and the function pblapply() which is a basically a
# parallel implementation of the mapping function
# lapply (or map from tidyverse) with a progress bar.

# When working on linux you can use forking and
# just map in parallel without any extra work
pblapply(grid_polygons$cell_id,
    function(cell) {
        grid_polygons %>%
            filter(cell_id == cell) %>%
            st_area()
    },
    cl = 2 # the number of cores
)

# When working on Windows you have to make a parallel 
# cluster first (number of cores in brackets)
cl <- makeCluster(2)

# then load any libraries you might need on the cluster
clusterEvalQ(cl, {
    library(sf)
    library(dplyr)
})

# and export any variables you need
clusterExport(cl, c("grid_polygons"))

# then run things in parallel 
pblapply(grid_polygons$cell_id,
    function(cell) {
        grid_polygons %>%
            filter(cell_id == cell) %>%
            st_area()
    },
    cl = cl # the name of your cluster
)

# once done you need to free up resouces by stopping 
# the cluster
stopCluster(cl)

# So basically if you can work on linux it's a lot easier!!

# Tip: when working with raster in parallel, then I often
# just use the filenames as the variable to iterate over
# I then load them in the parallel workers and write them
# out as files rather than returning them as objects
# That is more resource efficient and you don't need to 
# wrap the rast object from terra (see ?wrap)

# Hope this helps!

############
# Rectify/rotate data tryout
############
# rotated_data <- rast("./Aviris_data/strip_3328_cropped_rotated_2")
# plot()
# plot(ext(rotated_data))
# 
# plot(ext(sent_aoi_stack_crop))
# plot(ext(rotated_data), add=T)
# plotRGB(rotated_data, r=names(rotated_data)[54], g=names(rotated_data)[36], b=names(rotated_data)[20], stretch="lin", add=T)
# rectified_data <- rast("./Aviris_data/strip_3328_cropped_rect.envi")
# plot(ext(rectified_data), add=T)
# plotRGB(rectified_data, r=names(rectified_data)[54], g=names(rectified_data)[36], b=names(rectified_data)[20], stretch="lin", add=T)
# 
# rotated_data <- rast("./Aviris_data/strip_4159_cropped_2_rotated.dat")
# rectified_data <- rast("./Aviris_data/strip_4159_cropped_2_rect")
# # plot(ext(sent_aoi_stack_crop))
# plot(ext(rotated_data), add=T)
# plot(ext(rectified_data), add=T)
# 



