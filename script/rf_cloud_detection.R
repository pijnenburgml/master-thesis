library(rgdal)
library(terra)
library(caTools)
# install.packages("randomForest")
library(randomForest)
library(dplyr)
library(devtools)
library(gdalUtils)
library(sf)
library(exactextractr)


bg <- vect("~/scratch/data_rf_cloud_mask/strip_4159_background_annotation.shp")
plot(bg)
clouds <- vect("~/scratch/data_rf_cloud_mask/strip_4159_cloud_annotation.shp")
# clouds <- st_read("~/scratch/data_rf_cloud_mask/strip_4159_cloud_annotation.shp")
# clouds_path <- "~/scratch/data_rf_cloud_mask/strip_4159_cloud_annotation.shp"
plot(clouds, add=T, col="blue")

setwd("~/scratch/")
strip_4159_small_RF <- rast("~/scratch/strip_4159_small_RF_2_rotated.dat")
cloud_pix <- exact_extract(strip_4159_small_RF, clouds)
cloud_pix_terra <- terra::extract(strip_4159_small_RF, clouds)
bg_pix_terra <- terra::extract(strip_4159_small_RF, bg)
# write.table(cloud_pix_terra, file="~/data/output/strip_4159_cloud_pix")
# write.table(bg_pix_terra, file="~/data/output/strip_4159_bg_pix")
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
RF_whole_dataset <- cbind(cloud_pix_RF_t, bg_pix_RF_t)

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
pred_test <- predict(RF_model, newdata = validation_data, type= "class")
table(pred_test, validation_data$ID)  

# strip_4159_small_RF_path <- "~/scratch/data_rf_cloud_mask/stip_4159_small_RF.hdr"
# cloud_pix <- "~/scratch/data_rf_cloud_mask/cloud_pix_gdal.hdr"
# sf::gdal_utils(util="warp",
#                source= strip_4159_small_RF_path,
#                destination= cloud_pix,
#                options=c("cutline",  clouds_path, "of", "ENVI"))
# # GDAL Error 1: Cannot guess driver for ~/scratch/data_rf_cloud_mask/cloud_pix_gdal.ENVI
# 
# sf::gdal_utils(util="warp",
#                source= strip_4159_small_RF_path,
#                destination= cloud_pix,
#                options=c("cutline",  clouds_path))
# 
# # "crop to cutline"
# 
# 
# gdalwarp(srcfile = strip_4159_small_RF_path,
#                dstfile = cloud_pix, cutline = clouds_path)
# 
# gdal_cmd <- sprintf("gdalwarp -cutline %s", clouds_path, strip_4159_small_RF_path, cloud_pix)
# 
# 
# -cutline <datasource>

strip_4159_small_RF_rect <- rectify(strip_4159_small_RF)
writeRaster(strip_4159_small_RF_rect, filename = "strip_4159_small_RF_rect.envi")

strip_4159_half1 <- rast("~/scratch/data_rf_cloud_mask/strip_4159_half_1")
cloud_pix <- terra::extract(strip_4159_small_RF_rect, clouds)
strip_4159_half1 <- rectify(strip_4159_half1)
strip_4159_half2<- rast("strip_4159_half_2")
strip_4159_half2<- rectify(strip_4159_half2)
strip_4159 <- mosaic(strip_4159_half1, strip_4159_half2)
writeRaster(strip_4159_small_RF, filename ="strip_4159_small_RF_rect.envi")

cloud_pix <- extract(x=strip_4159_small_RF, y=clouds) %>% sample(size=10000, replace=F)
bg_pix <- extract(x=strip_4159, y=bg) %>% sample(size=10000, replace=F)

cloud_pix_half2 <- extract(x=)
bg_pix_half2 <- extract()

##########################
# Apply the random forest on another fligth strip 
##########################

strip_0747 <- rast("~/scratch/reflectance_data/ang20190801t160747_rfl")
strip_0747_cloud_class <- terra::predict(strip_0747, model=RF_model)









