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
################################################
# rectify envi file + save in the right format
################################################
non_rect <- rast("~/scratch/strip_4159_cropped_2")
yes_rect <- rectify(non_rect)
# writeRaster(yes_rect, filename = "~/data/Aviris_data/strip_4159_cropped_2_rect", filetype="ENVI", gdal=c("INTERLEAVE=BIL"), overwrite=T)

non_rect <- rast("~/scratch/strip_4159_cropped_2_merged_wa_sh_mask")
yes_rect <- rectify(non_rect)
writeRaster(yes_rect, filename = "~/data/Aviris_data/strip_4159_cropped_2_merged_wa_sh_mask_rect", filetype="ENVI", gdal=c("INTERLEAVE=BIL"), overwrite=T, datatype="INT1U")

################################################################################
# Variable definition
################################################################################
Datadir <- "Aviris_data"
# NameRaster <- "strip_3328_cropped_rect.envi"
NameRaster <- "strip_4159_cropped_2_rect"
# Define path for image file to be processed
Input_Image_File <- file.path(Datadir,NameRaster)
# Define path for corresponding mask file
# NameMask <- "strip_3328_cropped_1_merged_wa_sh_mask_rect.envi"
NameMask <- "strip_4159_cropped_2_merged_wa_sh_mask_rect"
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
nbCPU <- 4
MaxRAM <- 4
# number of clusters (spectral species)
nbclusters <- 50
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

# Select components from the PCA/SPCA/MNF raster
# Sel_PC = path of the file where selected components are stored
# Sel_PC <- select_PCA_components(Input_Image_File = Input_Image_File,
#                                 Output_Dir = Output_Dir,
#                                 PCA_Files = PCA_Output$PCA_Files,
#                                 TypePCA = PCA_Output$TypePCA,
#                                 File_Open = TRUE)
# 

SelectedPCs = c(1,2,7)

################################################################################
##                  Perform Spectral species mapping                          ##
## https://jbferet.github.io/biodivMapR/articles/biodivMapR_5.html            ##
################################################################################
print("MAP SPECTRAL SPECIES")
future.seed=TRUE
Kmeans_info <- map_spectral_species(Input_Image_File = Input_Image_File,
                                    Input_Mask_File = Input_Mask_File,
                                    Output_Dir = Output_Dir,
                                    SpectralSpace_Output = PCA_Output,
                                    nbclusters = nbclusters,
                                    nbCPU = nbCPU, MaxRAM = MaxRAM,
                                    SelectedPCs = SelectedPCs)


spectral_sp_map <- rast("~/data/biodivmapR_Aviris/RESULTS/strip_3328_cropped_rect/SPCA/SpectralSpecies/SpectralSpecies")
plot(spectral_sp_map$`Iter 2`, col=viridis_colors)

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
Index_Alpha   = c('Simpson')
Index_Alpha   = c('Shannon','Simpson')
map_alpha_div(Input_Image_File = Input_Image_File,
              Output_Dir = Output_Dir,
              TypePCA = TypePCA,
              window_size = window_size,
              nbCPU = nbCPU,
              MaxRAM = MaxRAM,
              Index_Alpha = Index_Alpha,
              nbclusters = nbclusters)

path <- paste("~/data/biodivmapR_Aviris/RESULTS/", NameRaster, "/SPCA/ALPHA/Shannon_10", sep="")
Shannon_map <- rast(path)
par(mfrow=c(1,2))
plot(Shannon_map, col=viridis_colors)
Simpson_map <- rast("~/data/biodivmapR_Aviris/RESULTS/strip_3328_cropped_rect/SPCA/ALPHA/Shannon_SD_10")
plot(Simpson_map, col=viridis_colors)
shannon_mean_map <- rast("~/data/biodivmapR_sent/RESULTS/sent_crop_envi/SPCA/ALPHA/Shannon_10_MeanFilter")


viridis_colors <- viridis::plasma(20)
na_col <- grey(0.8, alpha=0.5)
m <- ggplot() +
  geom_spatraster(data = Shannon_map, na.rm = TRUE, aes(fill=Shannon_10))+
  # theme_map()+
  theme_minimal()+
  scale_fill_gradientn(colours = viridis_colors, na.value = na_col) +
  labs(fill = substitute(paste("Shannon index", italic("H'"))))+
  plot_annotation(
    title = "Map of estimates of plant Shannon diversity index",
    theme = theme(plot.title = element_text(hjust=0.4, size = 16, face = "bold")))

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

path <- paste("~/data/Aviris_data/", NameRaster, sep="")
input <- rast(path)
plotRGB(input, r=names(input)[54], g=names(input)[36], b=names(input)[20], stretch="lin")
input_scaled <-input*2500

s <- ggplot() +
  geom_spatraster_rgb(data = input_scaled, r=58, g=37, b=23, interpolate=T)+
  # theme_map()+
  # theme_minimal()+
  plot_annotation(
    title = "True colour map",
    theme = theme(plot.title = element_text(hjust=0.4, size = 16, face = "bold")))
s

map_truecol <- s+
  ggspatial::annotation_scale(location="bl", pad_x=unit(0.5, "in"),
                              pad_y = unit(0.4, "in"), style="ticks", line_col="black", text_col="black")+
  
  ggspatial::annotation_north_arrow(location="bl", which_north=T, pad_x=unit(0.5, "in"),
                                    pad_y = unit(0.65, "in"), height = unit(0.6, "cm"), width=unit(0.6, "cm"))



cowplot::plot_grid(map_div, map_truecol, ncol = 2, align = "hv")



