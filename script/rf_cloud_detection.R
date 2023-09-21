library(rgdal)
library(terra)
library(caTools)
# install.packages("randomForest")
library(randomForest)
library(dplyr)
library(devtools)
library(sf)
library(ggplot2)

grid_polygons <- st_read("~/data/output/tilling/grid_polygons_aoi.shp")
grid_polygons


#############
# strip 5049
#############

# cl_5049 <- rast("~/scratch/data_rf_cloud_mask/strip_5049_clouds")
# cl_5049_rect <- rectify(cl_5049)
cl_5049 <- rast("~/scratch/data_rf_cloud_mask/strip_5049_clouds_rect.envi")
cl_5049_poly <- vect("~/scratch/data_rf_cloud_mask/strip_5049_cloud_polygon.shp")
cloud_pix_5049 <- terra::extract(cl_5049, cl_5049_poly)
write.table(cloud_pix_5049, file="~/scratch/data_rf_cloud_mask/strip_5049_cloud_pix")
# cloud_pix_5049 <- spatSample(cl_5049, size=10, replace=F, na.rm=T)
# to heavy to run
cloud_pix_5049 <- read.table("~/scratch/data_rf_cloud_mask/strip_5049_cloud_pix")
table(is.na(cloud_pix_5049))
set.seed(5)
cloud_pix_5049_idx <- sample(1:nrow(cloud_pix_5049), size=10000, replace=F)
cloud_pix_5049_sample <- cloud_pix_5049[cloud_pix_5049_idx, ]
cloud_pix_5049_t <- t(cloud_pix_5049_sample)
cloud_pix_5049_t["ID",] <- "cl"
rownames(cloud_pix_5049_t) <- c("ID", 1:425)
# write.table(cloud_pix_5049_t, file="~/data/output/random_forest_cloud_detection/cloud_pix_5049_t")
cloud_pix_5049_t <- read.table("~/data/output/random_forest_cloud_detection/cloud_pix_5049_t")


cl_5917 <- vect("~/scratch/data_rf_cloud_mask/strip_5917_cloud.shp")
# strip 5917 not in the aoi

##########
# strip 0708
##########
cl_0708 <- vect("~/scratch/data_rf_cloud_mask/strip_0708_clouds.shp")
cl_0708_ext <- vect(ext(cl_0708))
crs(cl_0708_ext) <- crs(cl_0708)
cl_0708_sf <- st_as_sf(cl_0708_ext)

# vizualisation 
ggplot(data=grid_polygons) +
  geom_sf(
    # aes(fill = cell_id, alpha=0.1)
  ) +
  geom_sf_label(aes(label=cell_id), cex=0.8)+
  geom_sf(data=cl_0708_sf, aes(fill="red", alpha=0.1))


cells <- c()
for (x in 1:nrow(grid_polygons)) {
  i <- st_intersects(grid_polygons$geometry[x], cl_0708_sf)
  i <- as.numeric(i)
  if(is.na(i)==F) {cells[x] <- grid_polygons$cell_id[x]}
}
cells <- cells[is.na(cells)==F]

paste0(paste("^",cells, sep=""), collapse = "|")

s_paths <- list.files(path = "~/scratch/reflectance_data", pattern = paste0(paste("^",cells, sep=""), collapse = "|"), full.names = T)
s_paths <- s_paths[-c(grep(s_paths, pattern=".hdr"))]
s_paths <- s_paths[-c(grep(s_paths, pattern=".aux.xml"))]
s_paths <- s_paths[-c(grep(s_paths, pattern="mask", ignore.case = T))]

sf::gdal_utils("buildvrt",
               source = s_paths, 
               destination = "/home/mpijne/scratch/strip_0708_cloud.vrt"
) 

sf::gdal_utils("translate",
               source="/home/mpijne/scratch/strip_0708_cloud.vrt",
               destination = "/home/mpijne/scratch/reflectance_data/strip_0708_cloud_data",
               options = c(
                 "-of", "ENVI",
                 "-co", "INTERLEAVE=BIL"
               )
)


cloud_0708_data <- rast("~/scratch/reflectance_data/strip_0708_cloud_data")
cloud_pix_0708 <- terra::extract(cloud_0708_data, cl_0708)
# write.table(cloud_pix_0708, file="~/data/output/random_forest_cloud_detection/strip_0708_cloud_pix")
cloud_pix_0708 <- read.table("~/data/output/random_forest_cloud_detection/strip_0708_cloud_pix")
table(is.na(cloud_pix_0708))
set.seed(5)
cloud_pix_0708_idx <- sample(1:nrow(cloud_pix_0708), size=10000, replace=F)
cloud_pix_0708_sample <- cloud_pix_0708[cloud_pix_0708_idx, ]
cloud_pix_0708_t <- t(cloud_pix_0708_sample)
cloud_pix_0708_t["ID",] <- "cl"
rownames(cloud_pix_0708_t) <- c("ID", 1:425)

# background pixel
bg_0708 <- vect("~/scratch/data_rf_cloud_mask/strip_0708_bg_polygon.shp")
plot(bg_0708)
bg_0708_1 <- bg_0708[1]
bg_0708_2 <- bg_0708[2]
bg_0708_3 <- bg_0708[3]

bg_0708_1_ext <- vect(ext(bg_0708_1))
crs(bg_0708_1_ext) <- crs(bg_0708_1)
bg_0708_1_sf <- st_as_sf(bg_0708_1_ext)

# vizualisation 
ggplot(data=grid_polygons) +
  geom_sf(
    # aes(fill = cell_id, alpha=0.1)
  ) +
  geom_sf_label(aes(label=cell_id), cex=0.8)+
  geom_sf(data=bg_0708_1_sf, aes(fill="red", alpha=0.1))

# not worth it :)

bg_0708_2_ext <- vect(ext(bg_0708_2))
crs(bg_0708_2_ext) <- crs(bg_0708_2)
bg_0708_2_sf <- st_as_sf(bg_0708_2_ext)

# vizualisation 
ggplot(data=grid_polygons) +
  geom_sf(
    # aes(fill = cell_id, alpha=0.1)
  ) +
  geom_sf_label(aes(label=cell_id), cex=0.8)+
  geom_sf(data=bg_0708_2_sf, aes(fill="red", alpha=0.1))


bg_0708_3_ext <- vect(ext(bg_0708_3))
crs(bg_0708_3_ext) <- crs(bg_0708_3)
bg_0708_3_sf <- st_as_sf(bg_0708_3_ext)

# vizualisation 
ggplot(data=grid_polygons) +
  geom_sf(
    # aes(fill = cell_id, alpha=0.1)
  ) +
  geom_sf_label(aes(label=cell_id), cex=0.8)+
  geom_sf(data=bg_0708_3_sf, aes(fill="red", alpha=0.1))

cells <- c()
for (x in 1:nrow(grid_polygons)) {
  i <- st_intersects(grid_polygons$geometry[x], bg_0708_2_sf$geometry)
  i <- as.numeric(i)
  if(is.na(i)==F) {cells[x] <- grid_polygons$cell_id[x]}
}
cells <- cells[is.na(cells)==F]
paste0(paste("^",cells, sep=""), collapse = "|")
s_paths <- list.files(path = "~/scratch/reflectance_data", pattern = paste0(paste("^",cells, sep=""), collapse = "|"), full.names = T)
s_paths <- s_paths[-c(grep(s_paths, pattern=".hdr"))]
s_paths <- s_paths[-c(grep(s_paths, pattern=".aux.xml"))]
sf::gdal_utils("buildvrt",
               source = s_paths, 
               destination = "/home/mpijne/scratch/strip_0708_2_bg.vrt"
) 

sf::gdal_utils("translate",
               source="/home/mpijne/scratch/strip_0708_2_bg.vrt",
               destination = "/home/mpijne/scratch/reflectance_data/strip_0708_bg_data_2",
               options = c(
                 "-of", "ENVI",
                 "-co", "INTERLEAVE=BIL"
               )
)

strip_0708_bg_data_2 <- rast("/home/mpijne/scratch/reflectance_data/strip_0708_bg_data_2")
strip_0708_bg_2_pix <- terra::extract(strip_0708_bg_data_2, bg_0708)
# write.table(strip_0708_bg_2_pix, file="~/data/output/random_forest_cloud_detection/strip_0708_bg_2_pix")
strip_0708_bg_2_pix <- read.table("~/data/output/random_forest_cloud_detection/strip_0708_bg_2_pix")
table(is.na(strip_0708_bg_2_pix))
na_idx <- which(is.na(strip_0708_bg_2_pix$strip_0708_bg_data_2_1))
strip_0708_bg_2_pix_noNA <- strip_0708_bg_2_pix[-na_idx,]
set.seed(5)
bg_2_0708_idx <- sample(1:nrow(strip_0708_bg_2_pix_noNA), size=10000, replace=F)
bg_2_0708_idx_sample <- strip_0708_bg_2_pix_noNA[bg_2_0708_idx, ]
bg_2_0708_idx_t <- t(bg_2_0708_idx_sample)
bg_2_0708_idx_t["ID",] <- "bg"
rownames(bg_2_0708_idx_t) <- c("ID", 1:425)


# polygon 3
cells <- c()
for (x in 1:nrow(grid_polygons)) {
  i <- st_intersects(grid_polygons$geometry[x], bg_0708_3_sf$geometry)
  i <- as.numeric(i)
  if(is.na(i)==F) {cells[x] <- grid_polygons$cell_id[x]}
}
cells <- cells[is.na(cells)==F]
paste0(paste("^",cells, sep=""), collapse = "|")
s_paths <- list.files(path = "~/scratch/reflectance_data", pattern = paste0(paste("^",cells, sep=""), collapse = "|"), full.names = T)
s_paths <- s_paths[-c(grep(s_paths, pattern=".hdr"))]
s_paths <- s_paths[-c(grep(s_paths, pattern=".aux.xml"))]
s_paths <- s_paths[-c(grep(s_paths, pattern="y_19"))] # not worth
sf::gdal_utils("buildvrt",
               source = s_paths, 
               destination = "/home/mpijne/scratch/strip_0708_3_bg.vrt"
) 

sf::gdal_utils("translate",
               source="/home/mpijne/scratch/strip_0708_3_bg.vrt",
               destination = "/home/mpijne/scratch/reflectance_data/strip_0708_bg_data_3",
               options = c(
                 "-of", "ENVI",
                 "-co", "INTERLEAVE=BIL"
               )
)

strip_0708_bg_data_3 <- rast("/home/mpijne/scratch/reflectance_data/strip_0708_bg_data_3")
strip_0708_bg_3_pix <- terra::extract(strip_0708_bg_data_3, bg_0708)
# write.table(strip_0708_bg_3_pix, file="~/data/output/random_forest_cloud_detection/strip_0708_bg_3_pix")
strip_0708_bg_3_pix <- read.table("~/data/output/random_forest_cloud_detection/strip_0708_bg_3_pix")
table(is.na(strip_0708_bg_3_pix))
na_idx <- which(is.na(strip_0708_bg_3_pix$strip_0708_bg_data_3_1))
strip_0708_bg_3_pix_noNA <- strip_0708_bg_3_pix[-na_idx,]
set.seed(5)
bg_3_0708_idx <- sample(1:nrow(strip_0708_bg_3_pix_noNA), size=10000, replace=F)
bg_3_0708_idx_sample <- strip_0708_bg_3_pix_noNA[bg_3_0708_idx, ]
bg_3_0708_idx_t <- t(bg_3_0708_idx_sample)
bg_3_0708_idx_t["ID",] <- "bg"
rownames(bg_3_0708_idx_t) <- c("ID", 1:425)


##############
# strip 5917
##############
# background pixel
bg_5917 <- vect("~/scratch/data_rf_cloud_mask/strip_5917_bg_polygon.shp")
plot(bg_5917)
bg_5917_ext <- vect(ext(bg_5917))
crs(bg_5917_ext) <- crs(bg_5917)
bg_5917_sf <- st_as_sf(bg_5917_ext)

# vizualisation 
ggplot(data=grid_polygons) +
  geom_sf(
    # aes(fill = cell_id, alpha=0.1)
  ) +
  geom_sf_label(aes(label=cell_id), cex=0.8)+
  geom_sf(data=bg_5917_sf, aes(fill="red", alpha=0.1))
# outside of aoi


##############
# strip 3328
##############
# background pixel
bg_3328 <- vect("~/scratch/data_rf_cloud_mask/strip_3328_bg_polygon.shp")
plot(bg_3328)
bg_3328_ext <- vect(ext(bg_3328))
crs(bg_3328_ext) <- crs(bg_3328)
bg_3328_sf <- st_as_sf(bg_3328_ext)

# vizualisation 
ggplot(data=grid_polygons) +
  geom_sf(
    # aes(fill = cell_id, alpha=0.1)
  ) +
  geom_sf_label(aes(label=cell_id), cex=0.8)+
  geom_sf(data=bg_3328_sf, aes(fill="red", alpha=0.1))

cells <- c()
for (x in 1:nrow(grid_polygons)) {
  i <- st_intersects(grid_polygons$geometry[x], bg_3328_sf$geometry)
  i <- as.numeric(i)
  if(is.na(i)==F) {cells[x] <- grid_polygons$cell_id[x]}
}
cells <- cells[is.na(cells)==F]
paste0(paste("^",cells, sep=""), collapse = "|")
s_paths <- list.files(path = "~/scratch/reflectance_data", pattern = paste0(paste("^",cells, sep=""), collapse = "|"), full.names = T)
s_paths <- s_paths[-c(grep(s_paths, pattern=".hdr"))]
s_paths <- s_paths[-c(grep(s_paths, pattern=".aux.xml"))]
sf::gdal_utils("buildvrt",
               source = s_paths, 
               destination = "/home/mpijne/scratch/strip_3328_bg.vrt"
) 

sf::gdal_utils("translate",
               source="/home/mpijne/scratch/strip_3328_bg.vrt",
               destination = "/home/mpijne/scratch/reflectance_data/strip_3328_bg_data",
               options = c(
                 "-of", "ENVI",
                 "-co", "INTERLEAVE=BIL"
               )
)

strip_3328_bg_data <- rast("/home/mpijne/scratch/reflectance_data/strip_3328_bg_data")
strip_3328_bg_pix <- terra::extract(strip_3328_bg_data, bg_3328)
# write.table(strip_3328_bg_pix, file="~/data/output/random_forest_cloud_detection/strip_3328_bg_pix")
strip_3328_bg_pix <- read.table("~/data/output/random_forest_cloud_detection/strip_3328_bg_pix")
table(is.na(strip_3328_bg_pix))
set.seed(5)
bg_3328_idx <- sample(1:nrow(strip_3328_bg_pix), size=10000, replace=F)
bg_3328_sample <- strip_3328_bg_pix[bg_3328_idx, ]
bg_3328_sample_t <- t(bg_3328_sample)
bg_3328_sample_t["ID",] <- "bg"
rownames(bg_3328_sample_t) <- c("ID", 1:425)


##############
# strip 4159
##############
bg <- vect("~/scratch/data_rf_cloud_mask/strip_4159_background_annotation.shp")
plot(bg)

clouds <- vect("~/scratch/data_rf_cloud_mask/strip_4159_cloud_annotation.shp")
plot(clouds, add=T, col="blue")

setwd("~/scratch/")
strip_4159_small_RF <- rast("~/scratch/strip_4159_small_RF_2_rotated.dat")
cloud_pix <- exact_extract(strip_4159_small_RF, clouds)
cloud_pix_terra <- terra::extract(strip_4159_small_RF, clouds)
bg_pix_terra <- terra::extract(strip_4159_small_RF, bg)
# write.table(cloud_pix_terra, file="~/data/output/strip_4159_cloud_pix")
# write.table(bg_pix_terra, file="~/data/output/strip_4159_bg_pix")
# cloud_pix_terra <- read.table("~/data/output/strip_4159_cloud_pix")
# bg_pix_terra <- read.table("~/data/output/strip_4159_bg_pix") #?

cloud_pix_terra <- read.table("~/data/output/strip_4159_cloud_pix")
bg_pix_terra <- read.table("~/data/output/strip_4159_bg_pix") #?
table(is.na(cloud_pix_terra))
na_idx <- which(is.na(bg_pix_terra$Rotate..Band.1.stip_4159_small_RF_2...376.719576.Nanometers.))
bg_pix_terra_no_na <- bg_pix_terra[-na_idx,]
table(is.na(bg_pix_terra_no_na))
set.seed(5)
cloud_pix_RF_idx <- sample(1:nrow(cloud_pix_terra), size=10000, replace=F)
cloud_pix_RF <- cloud_pix_terra[cloud_pix_RF_idx, ]
bg_pix_RF_idx <- sample(1:nrow(bg_pix_terra_no_na), size=10000, replace=F)
bg_pix_RF <- bg_pix_terra_no_na[bg_pix_RF_idx,]

cloud_pix_RF_t <- t(cloud_pix_RF)
cloud_pix_RF_t["ID",] <- "cl"
bg_pix_RF_t <- t(bg_pix_RF)
bg_pix_RF_t["ID",] <- "bg"

###############
# bringing data together
###############
RF_whole_dataset <- cbind(cloud_pix_RF_t, bg_pix_RF_t)
rownames(RF_whole_dataset) <- c("ID", 1:425)
RF_whole_dataset <- cbind(RF_whole_dataset, cloud_pix_5049_t,cloud_pix_0708_t,bg_2_0708_idx_t,bg_3_0708_idx_t,bg_3328_sample_t)
write.csv(RF_whole_dataset, file="~/data/output/random_forest_cloud_detection/RF_dataset")

split_idx <- sample.split(colnames(RF_whole_dataset), SplitRatio = 0.7)

training_data <- RF_whole_dataset[,split_idx]

validation_data <- RF_whole_dataset[,!split_idx]
validation_data <- as.data.frame(t(validation_data))
validation_data$ID <- as.factor(validation_data$ID)
validation_data <- validation_data %>%
  mutate_at(vars(-1), ~ as.numeric(.))
class(validation_data[,2])
colnames(validation_data)[2:426] <- 1:425
training_data <- as.data.frame(t(training_data))
training_data$ID <- as.factor(training_data$ID)
training_data <- training_data %>%
  mutate_at(vars(-1), ~ as.numeric(.))
class(training_data[,2])
colnames(training_data)[2:426] <- 1:425
RF_model <- randomForest(x=training_data[,2:426], y=training_data$ID)
saveRDS(RF_model, "~/data/output/random_forest_cloud_detection/cloud_classifier_RF.RData")
# save(RF_model,file = "~/data/output/random_forest_cloud_detection/cloud_classifier_RF.RData")
rm(RF_model)

RF_model <- readRDS("~/data/output/random_forest_cloud_detection/cloud_classifier_RF.RData")
# load("~/data/output/random_forest_cloud_detection/cloud_classifier_RF.RData")
pred_test <- predict(RF_model, newdata = validation_data, type= "class")
table(pred_test, validation_data$ID)  

# strip_4159_small_RF_rect <- rectify(strip_4159_small_RF)
# writeRaster(strip_4159_small_RF_rect, filename = "strip_4159_small_RF_rect.envi")

# strip_4159_half1 <- rast("~/scratch/data_rf_cloud_mask/strip_4159_half_1")
# cloud_pix <- terra::extract(strip_4159_small_RF_rect, clouds)
# strip_4159_half1 <- rectify(strip_4159_half1)
# strip_4159_half2<- rast("strip_4159_half_2")
# strip_4159_half2<- rectify(strip_4159_half2)
# strip_4159 <- mosaic(strip_4159_half1, strip_4159_half2)
# writeRaster(strip_4159_small_RF, filename ="strip_4159_small_RF_rect.envi")
# cloud_pix <- extract(x=strip_4159_small_RF, y=clouds) %>% sample(size=10000, replace=F)
# bg_pix <- extract(x=strip_4159, y=bg) %>% sample(size=10000, replace=F)

##########################
# Apply the random forest on another fligth strip 
##########################
# 
# tile_x10y14 <- rast("~/scratch/reflectance_data/x_10_y_14")
# names(tile_x10y14) <- 1:425
# x10y14_cloud_class <- terra::predict(tile_x10y14, model=RF_model)
# plot(x10y14_cloud_class)
# par(mfrow=c(1,2))
# plotRGB(tile_x10y14, r=54, g=36, b=20, stretch="lin")
# plot(x10y14_cloud_class)






