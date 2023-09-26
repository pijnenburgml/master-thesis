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

################################################################################
# Variable definition
################################################################################
Datadir <- "/scratch/mpijne/reflectance_data"
# NameRaster <- "strip_3328_cropped_rect.envi"
NameRaster <- "ang20190802t220708_rfl_rect"
# Define path for image file to be processed
Input_Image_File <- file.path(Datadir,NameRaster)
# Define path for corresponding mask file
# NameMask <- "strip_3328_cropped_1_merged_wa_sh_mask_rect.envi"
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
nbCPU <- 20
MaxRAM <- 8
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

# save(PCA_Output, file="~/data/biodivmapR_Aviris/RESULTS/ang20190802t220708_rfl_rect/SPCA/PCA/PCA_Output.Rdata")
get(load("~/data/biodivmapR_Aviris/RESULTS/ang20190802t220708_rfl_rect/SPCA/PCA/PCA_Output.Rdata"))


# path for the updated mask
Input_Mask_File <- PCA_Output$MaskPath

# way of visualizing PC variance explained
var_exp <- (PCA_Output$PCA_model$sdev^2/sum(PCA_Output$PCA_model$sdev^2))*100
barplot(var_exp, names.arg = colnames(PCA_Output$PCA_model$x))
cumsum(var_exp)
barplot(var_exp[1:10], names.arg=colnames(PCA_Output$PCA_model$x)[1:10])
barplot(var_exp[10:20], names.arg=colnames(PCA_Output$PCA_model$x)[10:20])
barplot(var_exp[50:425], names.arg=colnames(PCA_Output$PCA_model$x)[50:425])
viridis_colors <- viridis::inferno(60)
path <- paste("~/data/biodivmapR_Aviris/RESULTS/", NameRaster, "/SPCA/PCA/OutputPCA_30_PCs", sep="")
PCs <- rast(path)
plot(PCs[[1:8]], col=viridis_colors, legend=F)
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

# save(Kmeans_info, file="~/data/biodivmapR_Aviris/RESULTS_50_clusters/ang20190802t220708_rfl_rect/SPCA/SpectralSpecies/Kmeans_info_PC12.Rdata")
# Kmeans_info <- get(load("~/data/biodivmapR_Aviris/RESULTS_50_clusters/ang20190802t220708_rfl_rect/SPCA/SpectralSpecies/Kmeans_info_PC12.Rdata"))

library(RColorBrewer)
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
par(mfrow=c(1,2))
col=c("black", col_vector[2:n])
plot(spectral_sp_map$`Iter 1`, col=col)
plot(spectral_sp_map$`Iter 2`, col=col)

################################################################################
##                Perform alpha and beta diversity mapping                    ##
## https://jbferet.github.io/biodivMapR/articles/biodivMapR_6.html            ##
################################################################################
print("MAP ALPHA DIVERSITY")
Index_Alpha   = c('Shannon')
# Index_Alpha   = c('Shannon','Simpson')

spectral_sp <- rast("~/data/biodivmapR_Aviris/RESULTS/ang20190802t220708_rfl_rect/SPCA/SpectralSpecies/SpectralSpecies1289")
temp_rast <- rast(ext(spectral_sp), resolution=5)
spectral_sp_resample <- resample(spectral_sp, temp_rast)
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

##########
# plotting
##########

viridis_colors <- viridis::plasma(20)
na_col <- grey(0.8, alpha=0.5)
plot(Shannon_map, col=viridis_colors)

m <- ggplot() +
  geom_spatraster(data = Shannon_map, na.rm = TRUE, aes(fill=Shannon_20_PC1289))+
  # theme_map()+
  theme_minimal()+
  scale_fill_gradientn(colours = viridis_colors, na.value = na_col) +
  labs(fill = substitute(paste("Shannon index", italic("H'"))))
  # plot_annotation(
    # title = "Map of estimates of plant Shannon diversity index",
    # theme = theme(plot.title = element_text(hjust=0.4, size = 16, face = "bold")))
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

path_true_col <- Input_Image_File
input <- rast(path_true_col)
plotRGB(input, r=names(input)[54], g=names(input)[36], b=names(input)[20], stretch="lin")
input_scaled <-input*2500
t <- ggplot() +
  geom_spatraster_rgb(data = input_scaled, r=58, g=37, b=23, interpolate=T)+
  # theme_map()+
  # theme_minimal()+
  plot_annotation(
    title = "True colour map",
    theme = theme(plot.title = element_text(hjust=0.4, size = 16, face = "bold")))
t

map_truecol <- t+
  ggspatial::annotation_scale(location="bl", pad_x=unit(0.5, "in"),
                              pad_y = unit(0.4, "in"), style="ticks", line_col="black", text_col="black")+
  
  ggspatial::annotation_north_arrow(location="bl", which_north=T, pad_x=unit(0.5, "in"),
                                    pad_y = unit(0.65, "in"), height = unit(0.6, "cm"), width=unit(0.6, "cm"))



cowplot::plot_grid(map_div, map_truecol, ncol = 2, align = "hv")



