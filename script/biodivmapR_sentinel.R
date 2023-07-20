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
################################################################################
Datadir <- "~/data/biodivmapR_sent"
NameRaster <- "sent_crop_envi"
# Define path for image file to be processed
Input_Image_File <- file.path(Datadir,NameRaster)
# Define path for corresponding mask file
NameMask <- "mask_sent2_final_NA"
Input_Mask_File <- file.path(Datadir, NameMask)
# Input_Mask_File <- F
# Define path for master output directory where files produced during the process are saved
Output_Dir <- '~/data/biodivmapR_sent/RESULTS'
Output_Dir <- '~/data/biodivmapR_sent/RESULTS_cluster_20'
# dir.create(path = Output_Dir,recursive = T,showWarnings = F)
# NDVI_Thresh <- 0.8
# Blue_Thresh <- 500
# NIR_Thresh <- 1500
# Apply normalization with continuum removal?
Continuum_Removal <- F
# Type of dimensionality reduction
TypePCA <- 'SPCA'
# PCA FILTERING:        Set to TRUE if you want second filtering based on PCA outliers to be processed.
# Slower process
# Automatically set to FALSE if TypePCA     = 'MNF'
FilterPCA <- F
# window size forcomputation of spectral diversity
window_size <- 10
# # computational parameters
nbCPU <- 4
MaxRAM <- 4
# number of clusters (spectral species)
nbclusters <- 50
nbclusters <- 20

################################################################################
##                      Perform radiometric filtering                         ##
## https://jbferet.github.io/biodivMapR/articles/biodivMapR_3.html            ##
################################################################################
# print("PERFORM RADIOMETRIC FILTERING")
# Input_Mask_File <- perform_radiometric_filtering(Image_Path = Input_Image_File,
#                                                  Mask_Path = Input_Mask_File,
#                                                  Output_Dir = Output_Dir,
#                                                  TypePCA = TypePCA,
#                                                  NDVI_Thresh = NDVI_Thresh,
#                                                  Blue_Thresh = Blue_Thresh,
#                                                  NIR_Thresh = NIR_Thresh)
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
screeplot(PCA_Output$PCA_model)
# visual assesment of PCs
viridis_colors <- viridis::inferno(60)
PCs <- rast("biodivmapR_sent/RESULTS/sent_crop_envi/SPCA/PCA/OutputPCA_10_PCs")
PCs <- rast("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi/SPCA/PCA/OutputPCA_10_PCs")
plot(PCs, col=viridis_colors, legend=F)

viridis_colors <- viridis::plasma(20)
na_col <- grey(0.8, alpha=0.5)

pdf(file = "~/data/output/PCs_sentinel.pdf") 
for (x in 1:nlyr(PCs)) {
    print(ggplot() +
        geom_spatraster(data = PCs[[3]], na.rm = TRUE, aes(fill="PC 1"))+
        theme_minimal()+
        scale_fill_gradientn(colours = viridis_colors, na.value = na_col))
}
dev.off()

pdf(file = "~/data/output/PCs_sentinel.pdf") 
m <- ggplot() +
  geom_spatraster(data = PCs, na.rm = TRUE, aes(fill=`PC 10`))+
  # theme_map()+
  theme_minimal()+
  scale_fill_gradientn(colours = viridis_colors, na.value = na_col) 
  # labs(fill = substitute(italic("H'")))
m  
dev.off()



# Select components from the PCA/SPCA/MNF raster
# Sel_PC = path of the file where selected components are stored
Sel_PC <- select_PCA_components(Input_Image_File = Input_Image_File,
                                Output_Dir = Output_Dir,
                                PCA_Files = PCA_Output$PCA_Files,
                                TypePCA = PCA_Output$TypePCA,
                                File_Open = TRUE)


SelectedPCs = c(1,2,6,9)

################################################################################
##                  Perform Spectral species mapping                          ##
## https://jbferet.github.io/biodivMapR/articles/biodivMapR_5.html            ##
################################################################################
print("MAP SPECTRAL SPECIES")
future.seed=TRUE
Kmeans_info <- map_spectral_species(Input_Image_File = Input_Image_File,
                                    Input_Mask_File = PCA_Output$MaskPath,
                                    Output_Dir = Output_Dir,
                                    SpectralSpace_Output = PCA_Output,
                                    nbclusters = nbclusters,
                                    nbCPU = nbCPU, MaxRAM = MaxRAM,
                                    SelectedPCs = SelectedPCs)


spectral_sp_map <- rast("~/data/biodivmapR_sent/RESULTS/sent_crop_envi/SPCA/SpectralSpecies/SpectralSpecies")
spectral_sp_map <- rast("~/data/biodivmapR_sent/RESULTS_cluster_20//sent_crop_envi/SPCA/SpectralSpecies/SpectralSpecies")

plot(spectral_sp_map$`Iter 2`, col=viridis_colors)

# compute_spectral_species(PCA_Path = "biodivMapR_Example/03_RESULTS/S2A_T33NUD_20180104_Subset/SPCA/PCA/OutputPCA_8_PCs", 
#                          Spectral_Species_Path ="biodivMapR_Example/03_RESULTS/S2A_T33NUD_20180104_Subset/SPCA/SpectralSpecies/SpectralSpecies2", 
#                          Input_Mask_File = PCA_Output$MaskPath,
#                          PC_Select = c(1,2,5), Kmeans_info=Kmeans_info, nbCPU = 1)
# argument "Location_RW" is missing, with no default

################################################################################
##                Perform alpha and beta diversity mapping                    ##
## https://jbferet.github.io/biodivMapR/articles/biodivMapR_6.html            ##
################################################################################
SSD_Dir = c("~/data/biodivmapR_sent/RESULTS/sent_crop_envi/SPCA/alpha_metrics")
# dir.create(path = SSD_Dir,recursive = T,showWarnings = F)
Index_Alpha <- c("Shannon", "Simpson", "Fischer")

alpha_metric <- compute_alpha_metrics(
  Spectral_Species_Path = Kmeans_info$SpectralSpecies,
  SSD_Dir = SSD_Dir, 
  window_size = 10,
  Input_Mask_File = Input_Mask_File,
  nbclusters = 50,
  Index_Alpha = Index_Alpha,
  nbCPU = nbCPU,
  MaxRAM = MaxRAM
)
plot(alpha_metric$Shannon)

spectral_sp_distribution <- rast("~/data/biodivmapR_sent/RESULTS/sent_crop_envi/SPCA/alpha_metrics/SpectralSpecies_Distribution")
plot(spectral_sp_distribution$SpectralSpecies_Distribution_6)
shannon_map <- rast(alpha_metric$Shannon)
plot(shannon_map)
simpson_map <- rast(alpha_metric$Simpson)
fischer_map <- rast(alpha_metric$Fisher)
plot(simpson_map)
# why is it only 0?????
plot(fischer_map)

# other way to get alpha diversity map, also give only 0 for simpson diversity

print("MAP ALPHA DIVERSITY")
Index_Alpha   = c('Shannon','Simpson')
Index_Alpha   = c('Simpson')
map_alpha_div(Input_Image_File = Input_Image_File,
              Output_Dir = Output_Dir,
              TypePCA = TypePCA,
              window_size = window_size,
              nbCPU = nbCPU,
              MaxRAM = MaxRAM,
              Index_Alpha = Index_Alpha,
              nbclusters = nbclusters)

Shannon_map <- rast("~/data/biodivmapR_sent/RESULTS/sent_crop_envi/SPCA/ALPHA/Shannon_10")
Shannon_map <- rast("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi/SPCA/ALPHA/Shannon_10")
par(mfrow=c(1,2))
plot(Shannon_map)
plot(shannon_map)
Simpson_map <- rast("~/data/biodivmapR_sent/RESULTS/sent_crop_envi/SPCA/ALPHA/Simpson_10")
Simpson_map <- rast("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi/SPCA/ALPHA/Simpson_10")
plot(simpson_map)
plot(Simpson_map)
shannon_mean_map <- rast("~/data/biodivmapR_sent/RESULTS/sent_crop_envi/SPCA/ALPHA/Shannon_10_MeanFilter")

# spectral_sp_distribution <- rast("biodivMapR_Example/03_RESULTS/S2A_T33NUD_20180104_Subset/SPCA/SpectralSpecies/SpectralSpecies_Distribution")
# plot(spectral_sp_distribution)
#read_ENVI_header("biodivMapR_Example/03_RESULTS/S2A_T33NUD_20180104_Subset/SPCA/SpectralSpecies/SpectralSpecies_Distribution.hdr")


viridis_colors <- viridis::plasma(20)
na_col <- grey(0.8, alpha=0.5)
m <- ggplot() +
  geom_spatraster(data = Shannon_map, na.rm = TRUE, aes(fill=Shannon_10))+
  # theme_map()+
  theme_minimal()+
  scale_fill_gradientn(colours = viridis_colors, na.value = na_col) +
  labs(fill = substitute(paste("Shannon index", italic("H'"))))
m  

m+
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

sent_view <- rast("~/data/output/sent_crop_view.tif")
plotRGB(sent_view, r=3, g=2, b=1, scale=10000, stretch="lin", smooth=T)
sent_view_scale <-sent_view/7

s <- ggplot() +
  geom_spatraster_rgb(data = sent_view_scale, interpolate=T, r=3, g=2, b=1)+
  # theme_map()+
  theme_minimal()
s

s+
  ggspatial::annotation_scale(location="bl", pad_x=unit(0.5, "in"),
  pad_y = unit(0.4, "in"), style="ticks", line_col="white", text_col="white")+
  
  ggspatial::annotation_north_arrow(location="bl", which_north=T, pad_x=unit(0, "in"),
  pad_y = unit(0.4, "in"), height = unit(0.6, "cm"), width=unit(0.6, "cm"))










scale_fill_gradientn(colours = terrain.colors(10))


Ss <- rast("~/data/biodivmapR_sent/RESULTS_cluster_20//sent_crop_envi/SPCA/SpectralSpecies/SpectralSpecies")
plot(Ss$`Iter 4`, col=viridis_colors)


library(RColorBrewer)
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
par(mfrow=c(1,2))
col=c("black", col_vector[2:n])
plot(Ss$`Iter 4`, col=col)
plot(Ss$`Iter 5`, col=col)