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

################################################
# rectify envi file + save in the right format
################################################

non_rect <- rast("~/scratch/strip_4159_cropped_2")
yes_rect <- rectify(non_rect)
# writeRaster(yes_rect, filename = "~/data/Aviris_data/strip_4159_cropped_2_rect", filetype="ENVI", gdal=c("INTERLEAVE=BIL"), overwrite=T)

non_rect <- rast("~/scratch/strip_4159_cropped_2_merged_wa_sh_mask")
yes_rect <- rectify(non_rect)
# writeRaster(yes_rect, filename = "~/data/Aviris_data/strip_4159_cropped_2_merged_wa_sh_mask_rect", filetype="ENVI", gdal=c("INTERLEAVE=BIL"), overwrite=T)

# input <- rast("~/data/Aviris_data/strip_4159_cropped_2_rect.envi")
# table(is.na(as.vector(input$`461.869576 Nanometers`)))
# input_mask <- rast("~/data/Aviris_data/strip_4159_cropped_2_merged_wa_sh_mask_rect.envi")
# table(as.data.frame(input_mask)==0)
# 
# input <- rast("~/data/Aviris_data/strip_3328_cropped_rect.envi")
# table(is.na(as.vector(input$`376.719576 Nanometers`)))

################################################################################
Datadir <- "~/data/Aviris_data"
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
Output_Dir <- '~/data/biodivmapR_Aviris/RESULTS'
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
Input_Mask_File <- "PCA_Output$MaskPath"

# way of visualizing PC variance explained
var_exp <- (PCA_Output$PCA_model$sdev^2/sum(PCA_Output$PCA_model$sdev^2))*100
barplot(var_exp, names.arg = colnames(PCA_Output$PCA_model$x))
cumsum(var_exp)
barplot(var_exp[1:10], names.arg=colnames(PCA_Output$PCA_model$x)[1:10])
barplot(var_exp[10:20], names.arg=colnames(PCA_Output$PCA_model$x)[10:20])
barplot(var_exp[50:425], names.arg=colnames(PCA_Output$PCA_model$x)[50:425])
screeplot(PCA_Output$PCA_model)
viridis_colors <- viridis::inferno(60)
path <- paste("~/data/biodivmapR_Aviris/RESULTS/", NameRaster, "/SPCA/PCA/OutputPCA_30_PCs", sep="")
PCs <- rast(path)
plot(PCs[[1:8]], col=viridis_colors, legend=F)
pdf("PCstrip3328_to30.pdf")
# raster_pdf(filename = "PCstrip3328_to30.pdf")
plot(PCs[[1:10]], col=viridis_colors, legend=F)
dev.off
# not working to save in pdf... 

# visual assesment of PCs


# Select components from the PCA/SPCA/MNF raster
# Sel_PC = path of the file where selected components are stored
Sel_PC <- select_PCA_components(Input_Image_File = Input_Image_File,
                                Output_Dir = Output_Dir,
                                PCA_Files = PCA_Output$PCA_Files,
                                TypePCA = PCA_Output$TypePCA,
                                File_Open = TRUE)


SelectedPCs = c(1,2)

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

Shannon_map <- rast("~/data/biodivmapR_Aviris/RESULTS/strip_3328_cropped_rect/SPCA/ALPHA/Shannon_10")
par(mfrow=c(1,2))
plot(Shannon_map, col=viridis_colors)
Simpson_map <- rast("~/data/biodivmapR_Aviris/RESULTS/strip_3328_cropped_rect/SPCA/ALPHA/Shannon_SD_10")
plot(Simpson_map, col=viridis_colors)
shannon_mean_map <- rast("~/data/biodivmapR_sent/RESULTS/sent_crop_envi/SPCA/ALPHA/Shannon_10_MeanFilter")

# spectral_sp_distribution <- rast("biodivMapR_Example/03_RESULTS/S2A_T33NUD_20180104_Subset/SPCA/SpectralSpecies/SpectralSpecies_Distribution")
# plot(spectral_sp_distribution)
#read_ENVI_header("biodivMapR_Example/03_RESULTS/S2A_T33NUD_20180104_Subset/SPCA/SpectralSpecies/SpectralSpecies_Distribution.hdr")



viridis_colors <- viridis::inferno(20)
Ss <- rast("~/data/biodivmapR_sent/RESULTS_cluster_20//sent_crop_envi/SPCA/SpectralSpecies/SpectralSpecies")
plot(Ss$`Iter 4`, col=viridis_colors)


library(RColorBrewer)
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
par(mfrow=c(1,2))
col=c("black", col_vector[2:n])
plot(spectral_sp_map$`Iter 1`, col=col)
plot(spectral_sp_map$`Iter 2`, col=col)
