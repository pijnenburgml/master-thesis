# Example script for generating and processing tile geometries
# Jakob Assmann jakob.assmann@uzh.ch 27 June 2023

# Dependencies
library(sf)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(terra)
library(tidyterra)
setwd("~/data")
area_interest <- st_read("areas_of_interest.gpkg")
area_interest_proj <- st_transform(area_interest,32613)
sent_aoi_stack_crop <- rast("output/sent_crop_view.tif")
ext_aoi <- ext(area_interest_proj)
ext_sent <- ext(sent_aoi_stack_crop)

target_area <- data.frame(x = c(476480, 517470, 517470, 476480),
                          y = c(7657730, 7657730,  7687110, 7687110)) %>%
  st_as_sf(
    coords = c("x", "y"),
    crs = 32613 # UTM Zone 12 -> so coordiantes are in metres
  ) %>% summarize() %>%
  st_convex_hull()

plotRGB(sent_aoi_stack_crop, r=3, g=2, b=1, scale=10000, stretch="lin")
plot(target_area, add=T, col="red")

# # Generate an arbitary area of interest on Victoria Island
# target_area <- data.frame(x = c(468416, 468416 + 12600, 468416  + 12600, 468416),
#         y = c(7674690, 7674690, 7674690 + 12600, 7674690 + 12600 )) %>%
#             st_as_sf(
#                 coords = c("x", "y"),
#                 crs = 32613 # UTM Zone 12 -> so coordiantes are in metres
#             ) %>% summarize() %>%
#             st_convex_hull()

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
mask <- rast("~/master-thesis/sent_output/mask_sent2_final.tif")
ggplot() +
  geom_sf(
    data = grid_polygons,
    aes(fill = cell_id)
  ) +
  theme_map() +
  theme(legend.position = "none")
  # +geom_spatraster(data=mask, alpha=0.3)

##############################################################
boundary_dir <- "~/data/Aviris_data/boundaries"
boundaries <- list.files(path=boundary_dir, pattern="KML", full.names = T)
try <- st_read("~/data/Aviris_data/boundaries/ang20190801t143347_outline_KML.kml")
try_proj <- st_transform(try,32613)
st_intersects(grid_polygons[1,], try_proj$geometry)

ggplot(data=grid_polygons) +
  geom_sf(
    # aes(fill = cell_id, alpha=0.1)
  ) +
  geom_sf_label(aes(label=cell_id), cex=0.8)+
  geom_sf(data=try_proj, aes(fill="red"))
st_intersects(grid_polygons[grep(pattern="x_1_y_17", grid_polygons$cell_id),], try_proj)


# plotRGB(sent_aoi_stack_crop, r=3, g=2, b=1, scale=10000, stretch="lin")
# plot(try_proj$geometry, add=T, col="red")
# plot(grid_polygons$geometry, col=NA, border=1, )

boundaries_sf <- lapply(boundaries, FUN=st_read)
for(i in 1:length(boundaries_sf)){
  boundaries_sf[[i]] <- st_transform(boundaries_sf[[i]],32613) 
}

# substr(boundaries, 42, 75)
boundaries_sf[[1]]$Name

boundaries_name <- c()
for(i in 1:length(boundaries_sf)){
  boundaries_name[i] <- boundaries_sf[[i]]$Name
}
names(boundaries_sf) <- boundaries_name

try_multi <- c(boundaries_sf[[1]]$geometry, boundaries_sf[[2]]$geometry) %>% st_cast("MULTIPOLYGON")
names(try_multi) <- boundaries_name[1:2]
class(try_multi)
st_intersects(try_multi, grid_polygons[grep(pattern="x_5_y_17", grid_polygons$cell_id),], sparse = FALSE)
x <- st_intersects(try_multi, grid_polygons[grep(pattern="x_5_y_17", grid_polygons$cell_id),], sparse = T)
apply(x, 1, any)

p <- c()
for (y in 1:length(boundaries_sf)) {
  p[y] <- boundaries_sf[[y]]$geometry
}

cat(paste(rep("boundaries_sf[[", 16), 1:16, rep("]]$geometry", 16), sep=""), sep=",")
boundaries_multipolygon <- c(boundaries_sf[[1]]$geometry, boundaries_sf[[2]]$geometry,boundaries_sf[[3]]$geometry,boundaries_sf[[4]]$geometry,
  boundaries_sf[[5]]$geometry,boundaries_sf[[6]]$geometry,boundaries_sf[[7]]$geometry,boundaries_sf[[8]]$geometry,
  boundaries_sf[[9]]$geometry,boundaries_sf[[10]]$geometry,boundaries_sf[[11]]$geometry,boundaries_sf[[12]]$geometry,
  boundaries_sf[[13]]$geometry,boundaries_sf[[14]]$geometry,boundaries_sf[[15]]$geometry,boundaries_sf[[16]]$geometry)%>% st_cast("MULTIPOLYGON")
plot(boundaries_multipolygon)
names(boundaries_multipolygon) <- boundaries_name
idx <- which(st_intersects(boundaries_multipolygon,  grid_polygons[grep(pattern="x_5_y_17", grid_polygons$cell_id),], sparse = F)==T)
names(boundaries_multipolygon[1:2])

polygons_list <- list()
for (i in 1:length(boundaries_sf)) {
  polygons_list[[i]] <- boundaries_sf[[i]]$geometry[[1]]
}
polygons_list_XY <- st_zm(polygons_list)
names(polygons_list_XY) <- boundaries_name
boundaries_multipolygon <- st_multipolygon(x=polygons_list_XY, dim="XY")
plot(boundaries_multipolygon)

# polygons_list <- st_polygon()
# for (i in 1:length(boundaries_sf)) {
#   polygons_list[[i]] <- boundaries_sf[[i]]$geometry
# }
# boundaries_multipolygon <- st_cast(polygons_list, "MULTIPOLYGON")
# try_multi <- polygons_list %>% st_cast("MULTIPOLYGON")


geom <- st_geometry(grid_polygons[grep(pattern="x_5_y_17", grid_polygons$cell_id),])
id <- which(st_intersects(geom,boundaries_multipolygon))

for (x in 1:length(boundaries_sf)) {
  st_intersects(grid_polygons[grep(pattern="x_5_y_17", grid_polygons$cell_id),], boundaries_sf[[x]])  
}


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