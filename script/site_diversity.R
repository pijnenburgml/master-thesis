setwd("~/data/field_work/")
library(sf)
library(terra)
library(raster)
library(biodivMapR)

site_boundaries <-st_read("~/data/field_work/site_boundaries.gpkg")
plot_location <- st_read("~/data/field_work/plot_locations.gpkg")
presence_absence <- read.csv("~/data/field_work/presence_absence_cleaned.csv")
point_frame_data <- read.csv("~/data/field_work/pointframing.csv")
aoi <- st_read("~/data/field_work/areas_of_interest.gpkg")
aoi_proj <- st_transform(aoi,32613)
sent_view <- rast("~/data/output/sent_crop_view.tif")
plotRGB(sent_view, r=3, g=2, b=1, scale=1000, stretch="lin")
plot(site_boundaries, add=T, col="red")
plot(aoi_proj, add=T)
plot(aoi_proj)
plot(site_boundaries, add=T)
Shannon_map <- rast("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi/SPCA/ALPHA/Shannon_10")
plot(Shannon_map)
plot(site_boundaries, add=T)
site_boundaries[1,1] <- "Site 1"
plot(site_boundaries)
site_1 <- site_boundaries[1,]
site_2 <- site_boundaries[2,]
site_3 <- site_boundaries[3,]
plot(Shannon_map)
plot(site_1, add=T, col="red")

# st_write(site_1, "~/data/field_work/site_1.shp")
# st_write(site_2, "~/data/field_work/site_2.shp")
# st_write(site_3, "~/data/field_work/site_3.shp")

Sent_mask <- rast("~/data/biodivmapR_sent/mask_sent2_final_NA")

########
# biodivmapR from plot
########
# with 20 clusters

Datadir <- "~/data/biodivmapR_sent"
NameRaster <- "sent_crop_envi_BIL"
# Define path for image file to be processed
Input_Image_File <- file.path(Datadir,NameRaster)
# Define path for corresponding mask file
NameMask <- "mask_sent2_final_NA"
Input_Mask_File <- file.path(Datadir, NameMask)

# location of the directory where shapefiles used for validation are saved
VectorDir <- "~/data/field_work"
# list vector data
Path_Vector <- list.files("~/data/field_work", full.names = F, pattern=".shp")
Path_Vector <- list_shp(VectorDir)
Name_Vector <- tools::file_path_sans_ext(basename(Path_Vector))
Kmeans_info <- get(load("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/SpectralSpecies/Kmeans_info_PC1278.Rdata"))
PCA_Output <- get(load("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/PCA/PCA_Output.RData"))
# location of the spectral species raster needed for validation
Path_SpectralSpecies <- Kmeans_info$SpectralSpecies
nbclusters <- 20
# get diversity indicators corresponding to shapefiles (no partitioning of spectral dibversity based on field plots so far...)
Biodiv_Indicators <- diversity_from_plots_ML(Raster_SpectralSpecies = Path_SpectralSpecies, 
                                          Plots = Path_Vector,
                                          nbclusters = nbclusters)
# save(Biodiv_Indicators, file = "~/data/field_work/Biodive_Indicators_sent_20cluster_PC1278.Rdata")
Biodiv_Indicators <- get(load("~/data/field_work/Biodive_Indicators_sent_20cluster_PC1278.Rdata"))
Shannon_RS <- c(Biodiv_Indicators$Shannon)[[1]]

Spectral_sp_whole <- rast("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/SpectralSpecies/SpectralSpecies1278")
plot(Spectral_sp_whole$`Iter 1`)
#site 1
site_1 <- vect("~/data/field_work/site_1.shp")
Spectral_sp_whole_crop_s1 <- crop(Spectral_sp_whole, site_1)
mask_site_1 <- crop(Sent_mask, site_1)
plot(mask_site_1)
viridis_colors <- viridis::plasma(20)
terra::plot(Spectral_sp_whole_crop_s1$`Iter 1`, type="classes",col= viridis_colors)
table(as.data.frame(Spectral_sp_whole_crop_s1$`Iter 1`))
#site 2
site_2 <- vect("~/data/field_work/site_2.shp")
Spectral_sp_whole_crop_s2 <- crop(Spectral_sp_whole, site_2)
mask_site_2 <- crop(Sent_mask, site_2)
plot(mask_site_2)
viridis_colors <- viridis::plasma(20)
terra::plot(Spectral_sp_whole_crop_s2$`Iter 1`, type="classes",col= viridis_colors)
table(as.data.frame(Spectral_sp_whole_crop_s2$`Iter 1`))
#site 3
site_3 <- vect("~/data/field_work/site_3.shp")
Spectral_sp_whole_crop_s3 <- crop(Spectral_sp_whole, site_3)
mask_site_3 <- crop(Sent_mask, site_3)
plot(mask_site_3)
viridis_colors <- viridis::plasma(20)
terra::plot(Spectral_sp_whole_crop_s3$`Iter 1`, type="classes",col= viridis_colors)
table(as.data.frame(Spectral_sp_whole_crop_s3$`Iter 1`))


sent_view <- rast("~/data/output/sent_crop_view.tif")
sentinel_s1 <- crop(sent_view, site_1)
sentinel_s2 <- crop(sent_view, site_2)
sentinel_s3 <- crop(sent_view, site_3)

par(mfrow=c(2, 3))
terra::plotRGB(sentinel_s1, r=3, g=2, b=1, scale=1000, stretch="lin")
terra::plotRGB(sentinel_s2, r=3, g=2, b=1, scale=1000, stretch="lin")
terra::plotRGB(sentinel_s3, r=3, g=2, b=1, scale=1000, stretch="lin")
terra::plot(Spectral_sp_whole_crop_s1$`Iter 1`, type="classes",col= viridis_colors)
terra::plot(Spectral_sp_whole_crop_s2$`Iter 1`, type="classes",col= viridis_colors)
terra::plot(Spectral_sp_whole_crop_s3$`Iter 1`, type="classes",col= viridis_colors)

# #####
# # comparison between iteration
# #####
# # site 1
# pdf(file= "~/data/field_work/site_1_sp_sp_iterations_20clusters_PC129.pdf",   # The directory you want to save the file in
#     width = 15, # The width of the plot in inches
#     height = 15) # The height of the plot in inches
# 
# par(mfrow=c(3,4))
# 
# for (x in 1:nlyr(Spectral_sp_whole_crop_s1)) {
#   terra::plot(Spectral_sp_whole_crop_s1[[x]], type="classes", col=viridis_colors, main=paste("site 1 - iteration", x))
# }
# 
# dev.off()
# 
# # site 2
# pdf(file= "~/data/field_work/site_2_sp_sp_iterations_20clusters_PC129.pdf",   # The directory you want to save the file in
#     width = 15, # The width of the plot in inches
#     height = 15) # The height of the plot in inches
# 
# par(mfrow=c(3,4))
# 
# 
# for (x in 1:nlyr(Spectral_sp_whole_crop_s2)) {
#   terra::plot(Spectral_sp_whole_crop_s2[[x]], type="classes", col=viridis_colors, main=paste("site 2 - iteration", x))
# }
# 
# dev.off()
# 
# # site 3
# pdf(file= "~/data/field_work/site_3_sp_sp_iterations_20clusters_PC129.pdf",   # The directory you want to save the file in
#     width = 15, # The width of the plot in inches
#     height = 15) # The height of the plot in inches
# 
# par(mfrow=c(3,4))
# 
# 
# for (x in 1:nlyr(Spectral_sp_whole_crop_s3)) {
#   terra::plot(Spectral_sp_whole_crop_s3[[x]], type="classes", col=viridis_colors, main=paste("site 3 - iteration", x))
# }
# 
# dev.off()

library(ggplot2)
library(tidyterra)
library(cowplot)
sentinel_s1_scaled <- sentinel_s1/7
sentinel_s2_scaled <- sentinel_s2/7
sentinel_s3_scaled <- sentinel_s3/7

v1 <- ggplot() +
  geom_spatraster_rgb(data = sentinel_s1_scaled, r=3, g=2, b=1, interpolate=F)+
  geom_sf()+
  coord_sf(label_axes=c("----"))+
  theme_cowplot()

v2 <- ggplot() +
  geom_spatraster_rgb(data = sentinel_s2_scaled, r=3, g=2, b=1, interpolate=F)+
  geom_sf()+
  coord_sf(label_axes=c("----"))+
  theme_cowplot()


v3 <- ggplot() +
  geom_spatraster_rgb(data = sentinel_s3_scaled, r=3, g=2, b=1, interpolate=F)+
  geom_sf()+
  coord_sf(label_axes=c("----"))+
  theme_cowplot()


plot_grid(v1,v2,v3, align = "v", ncol=1, nrow=3)

viridis_colors <- viridis::plasma(20)
na_col <- grey(0.8, alpha=0.5)
d1 <- ggplot()+
  geom_spatraster(data=Spectral_sp_whole_crop_s1$`Iter 1`)+
  scale_fill_gradientn(colours = c(na_col,viridis_colors), n.breaks=21, guide = "legend")+
  theme_cowplot()+
  theme(legend.position = "none")+
  geom_sf()+
  coord_sf(label_axes=c("----"))

d2 <- ggplot()+
  geom_spatraster(data=Spectral_sp_whole_crop_s2$`Iter 1`)+
  theme_cowplot()+
  theme(legend.position = "none")+
  scale_fill_gradientn(colours = viridis_colors, na.value = na_col, n.breaks=21, guide = "legend")+
  geom_sf()+
  coord_sf(label_axes=c("----"))


d3 <- ggplot()+
  geom_spatraster(data=Spectral_sp_whole_crop_s3$`Iter 1`)+
  theme_cowplot()+
  # annotate("text", label = "Text No. 1")+
  scale_fill_gradientn(colours = viridis_colors, na.value = na_col, n.breaks=21, guide = "legend")+
  geom_sf()+
  coord_sf(label_axes=c("----"))+
  theme(legend.position = "none"
        # panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)
        # plot.margin=unit(c(0, 0, 0, 0), "cm")
        )
  

plot_grid(d1,d2,d3, align = "v", ncol=1, nrow=3)

plot_grid(v1, v2, v3, d1, d2, d3, ncol=2, nrow = 3,  byrow=F, labels=c("Site 1", "", "Site 2", "","Site 3"))

Biodiv_Indicators$Richness
Biodiv_Indicators$Shannon
