setwd("~/data")
library(sf)
library(tidyverse)
library(dplyr)
library(terra)
library(tidyterra)
library(stars)
library(raster)
library(remotes)
# remotes::install_github('cran/dissUtils')
library(dissUtils)
# remotes::install_github('jbferet/biodivMapR')
library(biodivMapR)
library(tools)
library(patchwork)
library(cowplot)

################################################################################
# Variable definition
################################################################################
Datadir <- "/scratch/mpijne/reflectance_data"
NameRaster <- "ang20190802t220708_rfl_rect"
# Define path for image file to be processed
Input_Image_File <- file.path(Datadir,NameRaster)
# Define path for corresponding mask file
NameMask <- "strip_0708_aoi_mask"
Input_Mask_File <- file.path(Datadir, NameMask)
# Input_Mask_File <- F
# Define path for master output directory where files produced during the process are saved
Output_Dir <- "biodivmapR_Aviris/RESULTS"
# dir.create(path = Output_Dir,recursive = T,showWarnings = F)
# Apply normalization with continuum removal?
Continuum_Removal <- F
# Type of dimensionality reduction
TypePCA <- 'SPCA'
# PCA FILTERING:        Set to TRUE if you want second filtering based on PCA outliers to be processed.
# Slower process
# Automatically set to FALSE if TypePCA     = 'MNF'
FilterPCA <- F
# window size for computation of spectral diversity
window_size <- 20
# # computational parameters
nbCPU <- 32
MaxRAM <- 8
# number of clusters (spectral species)
nbclusters <- 20

################################################################################
##                  Perform PCA & Dimensionality reduction                    ##
## https://jbferet.github.io/biodivMapR/articles/biodivMapR_4.html            ##
################################################################################
print("PERFORM DIMENSIONALITY REDUCTION")
PCA_Output <- perform_PCA(Input_Image_File = Input_Image_File,
                          Input_Mask_File = Input_Mask_File,
                          Output_Dir = Output_Dir,
                          TypePCA = TypePCA,
                          FilterPCA = FilterPCA,
                          Continuum_Removal = Continuum_Removal,
                          nbCPU = nbCPU,
                          MaxRAM = MaxRAM)
save(PCA_Output, file="~/data/biodivmapR_Aviris/RESULTS/ang20190802t220708_rfl_rect/SPCA/PCA/PCA_Output.Rdata", overwrite=T)
get(load("~/data/biodivmapR_Aviris/RESULTS/ang20190802t220708_rfl_rect/SPCA/PCA/PCA_Output.Rdata"))

# path for the updated mask
Input_Mask_File <- PCA_Output$MaskPath

# way of visualizing PC variance explained
var_exp <- (PCA_Output$PCA_model$sdev^2/sum(PCA_Output$PCA_model$sdev^2))*100
barplot(var_exp, names.arg = colnames(PCA_Output$PCA_model$x))
# visual assesment of PCs to be done in ArcGis
SelectedPCs = c(1,2,8,9)

################################################################################
##                  Perform Spectral species mapping                          ##
## https://jbferet.github.io/biodivMapR/articles/biodivMapR_5.html            ##
################################################################################
print("MAP SPECTRAL SPECIES")
future.seed=TRUE
Kmeans_info <- map_spectral_species_ML(Input_Image_File = Input_Image_File,
                                       Input_Mask_File = PCA_Output$MaskPath,
                                       Output_Dir = Output_Dir,
                                       SpectralSpace_Output = PCA_Output,
                                       nbclusters = nbclusters,
                                       nbCPU = nbCPU, MaxRAM = MaxRAM,
                                       SelectedPCs = SelectedPCs)
# save(Kmeans_info, file="~/data/biodivmapR_Aviris/RESULTS/ang20190802t220708_rfl_rect/SPCA/SpectralSpecies/Kmeans_info_PC1289.Rdata")
Kmeans_info <- get(load("~/data/biodivmapR_Aviris/RESULTS/ang20190802t220708_rfl_rect/SPCA/SpectralSpecies/Kmeans_info_PC1289.Rdata"))

################################################################################
##                Perform alpha and beta diversity mapping                    ##
## https://jbferet.github.io/biodivMapR/articles/biodivMapR_6.html            ##
################################################################################
print("MAP ALPHA DIVERSITY")
Index_Alpha   = c('Shannon')

# resample the spectral species map to exactly 5 m resolution
spectral_sp <- rast("~/data/biodivmapR_Aviris/RESULTS/ang20190802t220708_rfl_rect/SPCA/SpectralSpecies/SpectralSpecies1289")
temp_rast <- rast(ext(spectral_sp), resolution=5)
spectral_sp_resample <- resample(spectral_sp, temp_rast, method="near")
# writeRaster(spectral_sp_resample, filetype="ENVI", gdal="INTERLEAVE=BIL", overwrite=T, datatype="INT1U", filename="~/data/biodivmapR_Aviris/RESULTS/ang20190802t220708_rfl_rect/SPCA/SpectralSpecies/SpectralSpecies1289")

map_alpha_div_ML(Input_Image_File = Input_Image_File,
                 Output_Dir = Output_Dir,
                 TypePCA = TypePCA,
                 window_size = window_size,
                 nbCPU = nbCPU,
                 MaxRAM = MaxRAM,
                 Index_Alpha = Index_Alpha,
                 nbclusters = nbclusters, SelectedPCs = SelectedPCs)
path <- paste(Output_Dir, "/", NameRaster, "/SPCA/ALPHA/Shannon_", window_size, "_PC", paste0(SelectedPCs, collapse = ""), sep="")
Shannon_map <- rast(path)
plot(Shannon_map)

#########
# Correlation between sentinel and aviris Shannon index map of spectral species diversity
#########
library(spatialEco)
aviris_map <- rast("~/data/biodivmapR_Aviris/RESULTS/ang20190802t220708_rfl_rect/SPCA/ALPHA/Shannon_20_PC1289")
sentinel_map <- rast("~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_10_PC1278")
sentinel_map_crop <- crop(sentinel_map, aviris_map)
ext(aviris_map) <- ext(sentinel_map_crop)
corr_aviris_sent <- raster.modified.ttest(aviris_map, sentinel_map_crop)
viridis_colors <- viridis::plasma(20)
plot(corr_aviris_sent$corr, type="interval", breaks=c(-1,0,0.5,1), main="correlation between Shannon index based on aviris data and Sentinel-2 data", 
     col=c(viridis_colors[c(4,16)], "#009200"), plg=list(title="correlation coefficient"))
corr_aviris_sent$bin <- raster::cut(as.vector(corr_aviris_sent$corr), breaks = c(-1,0,0.5,1), labels=c("neg", "low", "high"))

corr_pix <- as.vector(corr_aviris_sent$bin) # assuming 1 --> neg, 2 --> low, 3 --> high
table(corr_pix)
corr_pix_noNA <- corr_pix[is.na(corr_pix)==F]
table(is.na(corr_pix_noNA))
perc_area_neg <- (table(corr_pix_noNA)[1]*10000)/(length(corr_pix_noNA)*10000)*100
perc_area_low <- (table(corr_pix_noNA)[2]*10000)/(length(corr_pix_noNA)*10000)*100
perc_area_high<- (table(corr_pix_noNA)[3]*10000)/(length(corr_pix_noNA)*10000)*100
perc_area_neg + perc_area_low


##########
# plotting
##########

# Shannon index derived from AVIRIS-NG data 

viridis_colors <- viridis::plasma(20)
na_col <- grey(0.8, alpha=0.5)
na_col <- "white"
plot(Shannon_map, col=viridis_colors)

m <- ggplot() +
  geom_spatraster(data = Shannon_map, na.rm = TRUE, aes(fill=Shannon_20_PC1289))+
  theme_map()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  scale_fill_gradientn(colours = viridis_colors, na.value = na_col, n.breaks=3) +
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  labs(fill = substitute(paste("Shannon index")))
m  
map_div <- m+
  ggspatial::annotation_north_arrow(
    location = "bl", which_north = "true",
    pad_x = unit(0.5, "in"), pad_y = unit(0.65, "in"),
    height = unit(0.7, "cm"),
    width = unit(0.7, "cm")
  )+
  ggspatial::annotation_scale(
    location = "bl", pad_x = unit(0.5, "in"),
    pad_y = unit(0.4, "in"), 
    style="ticks"
  )
map_div

save_plot(filename="~/data/output/final_plot/Shannon_map_aviris_res100_PC1289.png", map_div, base_height = 4, bg="white")


# True color map from the flight strip

strip_0708 <- rast("~/scratch/reflectance_data/ang20190802t220708_rfl_rect")
plotRGB(strip_0708, r=names(strip_0708)[54], g=names(strip_0708)[36], b=names(strip_0708)[20], stretch="lin")
strip_0708_RGB <- strip_0708[[c(56, 36, 20)]]
strip_0708_RGB_scaled <-strip_0708_RGB*1800
t <- ggplot() +
  geom_spatraster_rgb(data = strip_0708_RGB_scaled, r=1, g=2, b=3, interpolate=T)+
  theme_map()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))
t

map_truecol <- t+
  ggspatial::annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.5, "in"), pad_y = unit(0.65, "in"),
                                    height = unit(0.7, "cm"),width = unit(0.7, "cm"))+
  ggspatial::annotation_scale(location = "bl", pad_x = unit(0.5, "in"), pad_y = unit(0.4, "in"),style="ticks")
map_truecol
save_plot(filename="~/data/output/final_plot/RGB_view_avris_ext.png", map_truecol, base_height = 4, bg="white")


# Correlation between Aviris and Sentine-2 derived Shannon index diversity map

corr_plot <- ggplot() +
  geom_spatraster(data = corr_aviris_sent, na.rm = TRUE, aes(fill=bin))+
  theme_map()+
  scale_fill_manual(values=c(viridis_colors[c(4,16)], "#009200"), na.value="white", labels=c("-1 - 0", "0 - 0.5", "0.5 - 1", ""), name="Correlation coefficient", guide = guide_legend(reverse = TRUE))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  ggspatial::annotation_north_arrow(location = "bl", which_north = "true",pad_x = unit(0.5, "in"), pad_y = unit(0.65, "in"),height = unit(0.7, "cm"),width = unit(0.7, "cm"))+
  ggspatial::annotation_scale(location = "bl", pad_x = unit(0.5, "in"),pad_y = unit(0.4, "in"),style="ticks")
corr_plot
save_plot(corr_plot, filename = "~/data/output/final_plot/correlation_aviris_sent_shannon.png", base_height=4, bg = "white")


# sentinel correspondance
Shannon_map <- rast("~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_10_PC1278")
Shannon_map_crop <- crop(Shannon_map, aviris_map)
s <- ggplot() +
  geom_spatraster(data = Shannon_map_crop, na.rm = TRUE, aes(fill=Shannon_10_PC1278))+ #need to change the fill variable
  theme_map()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  scale_fill_gradientn(colours = viridis_colors, na.value = na_col,n.breaks=3) +
  labs(fill = substitute(paste("Shannon index ")))
s 

s_div <- s+
  ggspatial::annotation_north_arrow(
    location = "bl", which_north = "true",
    pad_x = unit(0.5, "in"), pad_y = unit(0.65, "in"),
    height = unit(0.7, "cm"),
    width = unit(0.7, "cm")
  )+
  ggspatial::annotation_scale(
    location = "bl", pad_x = unit(0.5, "in"),
    pad_y = unit(0.4, "in"), 
    style="ticks"
  )
s_div
save_plot(s_div, filename = "~/data/output/final_plot/Shannon_sentinel_aviris_aoi.png", base_height=4, bg = "white")


# Putting plots together

aviris_plot_together <- plot_grid(map_truecol, map_div, s_div, corr_plot, align = "v", ncol=1, axis = "r")
aviris_plot_together
save_plot(aviris_plot_together, filename = "~/data/output/final_plot/aviris_plot_together_n_RGB.png", base_height = 15, bg="white")



# # measure of distribution of the Shannon index
# aviris_map <- raster::raster("~/data/biodivmapR_Aviris/RESULTS/ang20190802t220708_rfl_rect/SPCA/ALPHA/Shannon_20_PC1289")
# try <- strucDiv(aviris_map, wsl=3, fun=contrast)
# aviris_map$bin <- raster::cut(as.vector(aviris_map$Shannon_20_PC1289), breaks = c(0,1.5,2,3), labels=c("low", "medium", "high"))
# plot(aviris_map$bin)
# aviris_shannon_bin <- as.vector(aviris_map$bin) # assuming 1 --> neg, 2 --> low, 3 --> high
# aviris_shannon_bin_noNA <- aviris_shannon_bin[is.na(aviris_shannon_bin)==F]
# perc_area_neg <- (table(aviris_shannon_bin_noNA)[1]*10000)/(length(aviris_shannon_bin_noNA)*10000)*100
# perc_area_low <- (table(aviris_shannon_bin_noNA)[2]*10000)/(length(aviris_shannon_bin_noNA)*10000)*100
# perc_area_high<- (table(aviris_shannon_bin_noNA)[3]*10000)/(length(aviris_shannon_bin_noNA)*10000)*100
# 
# sentinel_map <- raster::raster("~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_10_PC1278")
# try_sent <- strucDiv(sentinel_map, wsl=3, fun=contrast)
# sentinel_map$bin <- raster::cut(as.vector(sentinel_map$Shannon_10_PC1278), breaks = c(0,1.5,2,3), labels=c("low", "medium", "high"))
# plot(sentinel_map$bin)
# sentinel_shannon_bin <- as.vector(sentinel_map$bin) # assuming 1 --> neg, 2 --> low, 3 --> high
# sentinel_shannon_bin_noNA <- sentinel_shannon_bin[is.na(sentinel_shannon_bin)==F]
# perc_area_neg <- (table(sentinel_shannon_bin_noNA)[1]*10000)/(length(sentinel_shannon_bin_noNA)*10000)*100
# perc_area_low <- (table(sentinel_shannon_bin_noNA)[2]*10000)/(length(sentinel_shannon_bin_noNA)*10000)*100
# perc_area_high <- (table(sentinel_shannon_bin_noNA)[3]*10000)/(length(sentinel_shannon_bin_noNA)*10000)*100



