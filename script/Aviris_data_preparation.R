##############
# doing from the fligth strip directly
###############
path <- "/home/mpijne/scratch/reflectance_data/ang20190802t220708_rfl"
# ang20190802t220708_rfl
boundary_dir <- "~/data/Aviris_data/boundaries"
boundaries <- list.files(path=boundary_dir, pattern="KML", full.names = T)
strip_0708 <- st_read(boundaries[12])
strip_0708_n <- st_zm(strip_0708[1], drop=T, what="ZM")
strip_0708_df <- as.data.frame(strip_0708_n) 
strip_0708_polygon <- as(strip_0708_n, "Spatial")
st_write(strip_0708_n, "~/scratch/reflectance_data/boundary_strip_0708.shp")
boundary_strip_0708 <- st_read("~/scratch/reflectance_data/boundary_strip_0708.shp")
plot(boundary_strip_0708)
bbox <- st_bbox(boundary_strip_0708)
boundary_strip_0708_proj <- st_transform(boundary_strip_0708,32613)
boundary_strip_0708_crop <- st_crop(boundary_strip_0708_proj, area_interest_proj)
plot(boundary_strip_0708_proj)
plot(area_interest_proj, add=T, col="orange")
plot(boundary_strip_0708_crop, add=T, col="red")
plot(boundary_strip_0708_crop)
plot(boundary_strip_0708_proj, add=T, col="red")
bbox <- st_bbox(ext(boundary_strip_0708_crop))
plot(bbox)

# rectification of the flight strip 
sf::gdal_utils("warp",
               source = path, 
               destination = paste(path, "rect",sep="_"),
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

cell <- "ang20190802t220708_rfl_rect"
tile <- rast(file.path("/scratch/mpijne/reflectance_data", cell))
plot(ext(tile))
plotRGB(tile, add=T, r=54, g=36, b=20, stretch="lin")
RF_model <-readRDS("~/data/output/random_forest_cloud_detection/cloud_classifier_RF.RData")
build_area_mask <- rast("~/data/output/build_are_mask.tif")
sent_aoi_stack_crop <- rast("output/sent_crop_view.tif")

NIR_average <- mean(tile[[c(86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105)]])
hist(NIR_average, breaks=seq(terra::minmax(NIR_average)[1],terra::minmax(NIR_average)[2]+0.05, by=0.01))
shade_mask <- ifel(NIR_average<0.05, 0, 1)
plot(shade_mask)
green_average <- mean(tile[[c(26, 27, 28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45)]])
NDWI <- (green_average-NIR_average)/(green_average+NIR_average)
hist(NDWI)
water_mask <- ifel(NDWI>0.1, 0, 1)
plot(water_mask)
names(tile) <- 1:425
library(randomForest)
cloud_class <- terra::predict(object=tile, model=RF_model, na.rm=T)
cloud_mask <- ifel(cloud_class=="cl", 0, 1)
plot(cloud_mask)
build_area_mask <- rast("~/data/output/build_are_mask.tif")
temp_rast <- rast(ext(tile), resolution = res(tile))
build_area_mask_5res <- resample(build_area_mask, temp_rast, method="near")
building_mask <- crop(build_area_mask_5res, tile)
mask <- mosaic(shade_mask, water_mask, cloud_mask, building_mask, fun="min")
plot(mask)
mask <- ifel(mask==0, NA, 1) #is it necessary?
# writeRaster(mask, filename = paste("/scratch/mpijne/reflectance_data/strip_0708_aoi_mask", sep=""),filetype="ENVI", gdal="INTERLEAVE=BSQ", overwrite=T, datatype="INT1U")
mask <- rast("/scratch/mpijne/reflectance_data/strip_0708_aoi_mask")
