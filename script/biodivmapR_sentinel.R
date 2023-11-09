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
library(ggplot2)
library(cowplot)
library(tools)
library(patchwork)
################################################################################
# save the data in the proper format
sent <- rast("~/data/biodivmapR_sent/sent_crop_envi")
writeRaster(sent, filename = "~/data/biodivmapR_sent/sent_crop_envi_BIL", filetype="ENVI", gdal=c("INTERLEAVE=BIL"), overwrite=T)

################################################################################

##########
# creating sentinel file corresponding to aviris flight strip 0708
##########
sent_whole <- rast("~/data/biodivmapR_sent/sent_crop_envi_BIL")
mask_whole <- rast("~/data/biodivmapR_sent/mask_sent2_final_NA")
mask_aviris <- rast("/scratch/mpijne/reflectance_data/strip_0708_aoi_mask")
plot(mask_aviris)
ext(mask_aviris)
plot(mask_whole)
ext(mask_whole)
temp_rast <- rast(ext(mask_whole), resolution=10)
mask_aviris_10 <- resample(mask_aviris, temp_rast, method="near")
ext(mask_aviris_10)
plot(mask_whole)
plot(mask_aviris_10, add=T)
sent_crop_mask <- mask(mask_whole, mask_aviris_10)
plot(sent_crop_mask, colNA="red")
# writeRaster(sent_crop_mask, filetype="ENVI", gdal="INTERLEAVE=BSQ", overwrite=T, datatype="INT1U", filename = "~/data/biodivmapR_sent/mask_matchin_aviris")


###########
# Setting parameters
###########

Datadir <- "~/data/biodivmapR_sent"
NameRaster <- "sent_crop_envi_BIL"
# Define path for image file to be processed
Input_Image_File <- file.path(Datadir,NameRaster)
# Define path for corresponding mask file
# => uncomment the right mask depending on the analysis
# NameMask <- "mask_sent2_final_NA"
# NameMask <- "mask_matchin_aviris"
Input_Mask_File <- file.path(Datadir, NameMask)
# Input_Mask_File <- F
# Define path for master output directory where files produced during the process are saved
# => un-comment the right directory
# Output_Dir <- '~/data/biodivmapR_sent/RESULTS_cluster_20'
# Output_Dir <- "~/data/biodivmapR_sent/RESULTS_PC_selection"
# Output_Dir <- "~/data/biodivmapR_sent/RESULTS_aviris_extend"
# dir.create(path = Output_Dir,recursive = T,showWarnings = F)
# Apply normalization with continuum removal?
Continuum_Removal <- F
# Type of dimensionality reduction
TypePCA <- 'SPCA'
# PCA FILTERING:        Set to TRUE if you want second filtering based on PCA outliers to be processed.
# Slower process
# Automatically set to FALSE if TypePCA     = 'MNF'
FilterPCA <- F
# window size forcomputation of spectral diversity
window_size <-10
# # computational parameters
nbCPU <- 6
MaxRAM <- 4
# number of clusters (spectral species)
nbclusters <- 50
nbclusters <- 20

################################################################################
##                  Perform PCA & Dimensionality reduction                    ##
## https://jbferet.github.io/biodivMapR/articles/biodivMapR_4.html            ##
################################################################################
print("PERFORM DIMENSIONALITY REDUCTION")
PCA_Output <- perform_PCA_ML(Input_Image_File = Input_Image_File,
                          Input_Mask_File = Input_Mask_File,
                          Output_Dir = Output_Dir,
                          TypePCA = TypePCA,
                          FilterPCA = FilterPCA,
                          Continuum_Removal = Continuum_Removal,
                          nbCPU = nbCPU,
                          MaxRAM = MaxRAM)
# save(PCA_Output, file="~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/PCA/PCA_Output.RData")
PCA_Output <- get(load("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/PCA/PCA_Output.RData"))

# save(PCA_Output, file="~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/PCA/PCA_Output.Rdata")
PCA_Output <- get(load("~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/PCA/PCA_Output.Rdata"))


# path for the updated mask
Input_Mask_File <- PCA_Output$MaskPath

# visualizing PC variance explained
screeplot(PCA_Output$PCA_model)

# Select PC based on visual assessment in ArcGIS software
SelectedPCs = c(1,2,7,8)

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
# save(Kmeans_info, file="~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/SpectralSpecies/Kmeans_info_PC1278.Rdata")
Kmeans_info <- get(load("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/SpectralSpecies/Kmeans_info_PC1278.Rdata"))

# save(Kmeans_info, file="~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/SpectralSpecies/Kmeans_info_PC1278.Rdata")
Kmeans_info <- get(load("~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/SpectralSpecies/Kmeans_info_PC1278.Rdata"))

# spectral_sp_map <- rast("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi/SPCA/SpectralSpecies/SpectralSpecies")
# spectral_sp_map <- rast("~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/SpectralSpecies/SpectralSpecies1278")
# spectral_sp_map
# viridis_colors <- viridis::plasma(20)
# plot(spectral_sp_map$`Iter 2`, col=viridis_colors)


################################################################################
##                Perform alpha and beta diversity mapping                    ##
## https://jbferet.github.io/biodivMapR/articles/biodivMapR_6.html            ##
################################################################################
print("MAP ALPHA DIVERSITY")
Index_Alpha   = c('Shannon')
map_alpha_div_ML(Input_Image_File = Input_Image_File,
              Output_Dir = Output_Dir,
              TypePCA = TypePCA,
              window_size = window_size,
              nbCPU = nbCPU,
              MaxRAM = MaxRAM,
              Index_Alpha = Index_Alpha,
              nbclusters = nbclusters, SelectedPCs = SelectedPCs)
path <- paste(Output_Dir, "/", NameRaster, "/SPCA/ALPHA/Shannon_", window_size, "_PC", paste0(SelectedPCs, collapse = ""),sep="")
Shannon_map <- rast(path)
# Shannon_map <- rast("~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_10_PC1278")
# Shannon_map <- rast("~/data/biodivmapR_sent/RESULTS/sent_crop_envi/SPCA/ALPHA/Shannon_10")
# Shannon_map <- rast("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_10")

# Plotting

viridis_colors <- viridis::plasma(20)
na_col <- grey(0.8, alpha=0.5)
na_col <- "white"
m <- ggplot() +
  geom_spatraster(data = Shannon_map, na.rm = TRUE, aes(fill=Shannon_10_PC1278))+ #need to change the fill variable
  theme_map()+
  # theme_minimal()+
  scale_fill_gradientn(colours = viridis_colors, na.value = na_col, n.breaks=3) +
  labs(fill = substitute(paste("Shannon index ")))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))
m 

map_div <- m+
  ggspatial::annotation_north_arrow(
    location = "bl", which_north = "true",
    pad_x = unit(0.5, "in"), pad_y = unit(0.45, "in"),
    height = unit(0.7, "cm"),
    width = unit(0.7, "cm")
  )+
  ggspatial::annotation_scale(
    location = "bl", pad_x = unit(0.5, "in"),
    pad_y = unit(0.2, "in"), 
    style="ticks"
  )
map_div
save_plot(filename = "~/data/output/final_plot/Shannon_map_sent_PC1278.png", base_height=4, map_div,bg = "white")
save_plot(filename = "~/data/output/final_plot/Shannon_map_sent_aviris_aoi_PC1278.svg", base_height=3, map_div,bg = "white")

# ######################
# # Beta diversity map
# ######################
# print("MAP BETA DIVERSITY")
# map_beta_div(Input_Image_File = Input_Image_File, 
#              Output_Dir = Output_Dir, 
#              TypePCA = TypePCA,
#              window_size = window_size, 
#              nbCPU = nbCPU, 
#              MaxRAM = MaxRAM,
#              nbclusters = nbclusters)
# 
# window_size <- 30
# path <- paste(Output_Dir, "/", NameRaster, "/SPCA/BETA/BetaDiversity_BCdiss_PCO_", window_size, sep="")
# Beta_map <- rast(path)
# Beta_map_10 <- rast(path)
# plotRGB(Beta_map, r=1, g=2, b=3, stretch="lin")
# 
# b <- ggplot() +
#   geom_spatraster_rgb(data = Beta_map, r=1, g=2, b=3, interpolate = T)+
#   # theme_map()+
#   theme_minimal()
#   # scale_fill_gradientn(colours = viridis_colors, na.value = na_col) +
#   # labs(fill = substitute(paste("Shannon index", italic("H'"))))
# b  


# #######################################
# # Comparison between Shannon index map produced with different combination of PC
# #######################################
# Shannon_PC12 <- rast("~/data/biodivmapR_sent/RESULTS_PC_selection/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_10_PC12")
# Shannon_PC129 <- rast("~/data/biodivmapR_sent/RESULTS_PC_selection/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_10_PC129")
# Shannon_PC12789 <- rast("~/data/biodivmapR_sent/RESULTS_PC_selection/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_10_PC12789")
# par(mfrow=c(1,2))
# plot(Shannon_PC12)
# plot(Shannon_PC129)
# plot(Shannon_PC12789)
# par(mfrow=c(1,1))
# 
# viridis_colors <- viridis::plasma(20)
# na_col <- grey(0.8, alpha=0.5)
# map_shannon_PC12 <- ggplot() +
#   geom_spatraster(data = Shannon_PC12, na.rm = TRUE, aes(fill=Shannon_10_PC12))+ #need to change the fill variable
#   # theme_map()+
#   theme_minimal()+
#   scale_fill_gradientn(colours = viridis_colors, na.value = na_col) +
#   labs(fill = substitute(paste("Shannon index", italic("H'"))))+
#   plot_annotation(
#     title = "Estimates of plant Shannon diversity index using PCs 1, 2",
#     theme = theme(plot.title = element_text(hjust=0.4, size = 16, face = "bold")))+
#   ggspatial::annotation_scale(location="bl", pad_x=unit(0.7, "in"),
#                               pad_y = unit(0.6, "in"), style="ticks", line_col="black", text_col="black")+
#   ggspatial::annotation_north_arrow(location="bl", which_north=T, pad_x=unit(0.7, "in"),
#                                     pad_y = unit(0.8, "in"), height = unit(0.6, "cm"), width=unit(0.6, "cm"))
# map_shannon_PC12  
# 
# 
# map_shannon_PC129 <- ggplot() +
#   geom_spatraster(data = Shannon_PC129, na.rm = TRUE, aes(fill=Shannon_10_PC129))+ #need to change the fill variable
#   # theme_map()+
#   theme_minimal()+
#   scale_fill_gradientn(colours = viridis_colors, na.value = na_col) +
#   labs(fill = substitute(paste("Shannon index", italic("H'"))))+
#   plot_annotation(
#     title = "Estimates of plant Shannon diversity index using PCs 1,2,9",
#     theme = theme(plot.title = element_text(hjust=0.4, size = 16, face = "bold")))+
#   ggspatial::annotation_scale(location="bl", pad_x=unit(0.7, "in"),
#                               pad_y = unit(0.6, "in"), style="ticks", line_col="black", text_col="black")+
#   ggspatial::annotation_north_arrow(location="bl", which_north=T, pad_x=unit(0.7, "in"),
#                                     pad_y = unit(0.8, "in"), height = unit(0.6, "cm"), width=unit(0.6, "cm"))
# map_shannon_PC129
# 
# cowplot::plot_grid(map_shannon_PC12, map_shannon_PC129, ncol = 2, align = "hv")
# library(spatialEco)
# library(SpatialPack)
# corr_diff_PC <- raster.modified.ttest(Shannon_PC12, Shannon_PC129)
# plot(corr_diff_PC$corr, type="interval", breakby="cases")
# writeRaster(corr_diff_PC, "~/data/biodivmapR_sent/RESULTS_PC_selection/sent_crop_envi_BIL/SPCA/ALPHA/correlation_PC12_PC129.envi")
# corr_diff_PC <- rast("~/data/biodivmapR_sent/RESULTS_PC_selection/sent_crop_envi_BIL/SPCA/ALPHA/correlation_PC12_PC129.envi")
# terra::plot(corr_diff_PC$corr, col=viridis_colors, type="continuous")
# plot(corr_diff_PC$corr, type="interval", breakby="cases", main="correlation between Shannon index produced with PC 1,2 and PC 1,2,9")

# map_shannon_PC12789 <- ggplot() +
#   geom_spatraster(data = Shannon_PC12789, na.rm = TRUE, aes(fill=Shannon_10_PC12789))+ #need to change the fill variable
#   # theme_map()+
#   theme_minimal()+
#   scale_fill_gradientn(colours = viridis_colors, na.value = na_col) +
#   labs(fill = substitute(paste("Shannon index", italic("H'"))))+
#   plot_annotation(
#     title = "Estimates of plant Shannon diversity index using PCs 1,2,7,8,9",
#     theme = theme(plot.title = element_text(hjust=0.4, size = 16, face = "bold")))+
#   ggspatial::annotation_scale(location="bl", pad_x=unit(0.7, "in"),
#                               pad_y = unit(0.6, "in"), style="ticks", line_col="black", text_col="black")+
#   ggspatial::annotation_north_arrow(location="bl", which_north=T, pad_x=unit(0.7, "in"),
#                                     pad_y = unit(0.8, "in"), height = unit(0.6, "cm"), width=unit(0.6, "cm"))
# map_shannon_PC12789  
# cowplot::plot_grid(map_shannon_PC12, map_shannon_PC129, map_shannon_PC12789, ncol = 2, align = "hv")
# cowplot::plot_grid(map_shannon_PC12, map_shannon_PC12789, ncol = 2, align = "hv")
# corr_diff_PC_12_12789 <- raster.modified.ttest(Shannon_PC12, Shannon_PC12789)
# # writeRaster(corr_diff_PC_12_12789, "~/data/biodivmapR_sent/RESULTS_PC_selection/sent_crop_envi_BIL/SPCA/ALPHA/correlation_PC12_PC12789.envi")
# corr_diff_PC_12_12789 <- rast("~/data/biodivmapR_sent/RESULTS_PC_selection/sent_crop_envi_BIL/SPCA/ALPHA/correlation_PC12_PC12789.envi")
# plot(corr_diff_PC_12_12789$corr, type="interval", breakby="eqint", main="correlation between Shannon index produced with PC 1,2 and PC 1,2,7,8,9", col=viridis_colors)
# plot(corr_diff_PC_12_12789$corr, type="interval", breaks=c(-1,-0.5,-0.1,0.1,0.5, 1), main="correlation between Shannon index produced with PC 1,2 and PC 1,2,7,8,9", col=viridis_colors[c(1,5,10,15,18)])
# plot(corr_diff_PC_12_12789$corr, type="interval", breaks=c(-1,-0.5,-0.1,0.1,0.5, 1), main="correlation between Shannon index produced with PC 1,2 and PC 1,2,7,8,9", col=c(viridis_colors[c(4,8,12,16)], "#009200"), plg=list(title="correlation coefficient"))

