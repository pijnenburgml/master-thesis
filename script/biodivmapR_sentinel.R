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
# sent <- rast("~/data/biodivmapR_sent/sent_crop_envi")
# writeRaster(sent, filename = "~/data/biodivmapR_sent/sent_crop_envi_BIL", filetype="ENVI", gdal=c("INTERLEAVE=BIL"), overwrite=T)

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
mask_aviris_10 <- resample(mask_aviris, temp_rast)
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
# NameMask <- "mask_sent2_final_NA"
NameMask <- "mask_matchin_aviris"
Input_Mask_File <- file.path(Datadir, NameMask)
# Input_Mask_File <- F
# Define path for master output directory where files produced during the process are saved
Output_Dir <- '~/data/biodivmapR_sent/RESULTS_cluster_50'
Output_Dir <- '~/data/biodivmapR_sent/RESULTS_cluster_20'
Output_Dir <- "~/data/biodivmapR_sent/RESULTS_PC_selection"
Output_Dir <- "~/data/biodivmapR_sent/RESULTS_aviris_extend"
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
window_size <-10
# # computational parameters
nbCPU <- 2
MaxRAM <- 6
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
# save(PCA_Output, file="~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/PCA/PCA_Output.RData")
PCA_Output <- get(load("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/PCA/PCA_Output.RData"))

# save(PCA_Output, file="~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/PCA/PCA_Output.Rdata")
PCA_Output <- get(load("~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/PCA/PCA_Output.Rdata"))


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


SelectedPCs = c(1,2,7,8)

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
# save(Kmeans_info, file="~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/SpectralSpecies/Kmeans_info_PC12.Rdata")
Kmeans_info <- get(load("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/SpectralSpecies/Kmeans_info.Rdata"))


# To have the name of the PCs selected to map the spectral species in the filename 
# run the function saved in the document map_spectral_sp_PC_naming.R to have in the 
# global environment

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
Kmeans_info <- get(load("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/SpectralSpecies/Kmeans_info_PC1278.Rdata"))


spectral_sp_map <- rast("~/data/biodivmapR_sent/RESULTS/sent_crop_envi/SPCA/SpectralSpecies/SpectralSpecies")
spectral_sp_map <- rast("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi/SPCA/SpectralSpecies/SpectralSpecies")

plot(spectral_sp_map$`Iter 2`, col=viridis_colors)

PCs <- paste0(SelectedPCs, collapse = "")
path <- file.path(Output_Dir, NameRaster, TypePCA, "SpectralSpecies", paste0("SpectralSpecies", PCs, collapse = ""))
spectral_sp <- rast(path)
plot(spectral_sp$`Iter 1`)
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

path <- file.path(Output_Dir, NameRaster, TypePCA, "SpectralSpecies/SpectralSpecies_Distribution")
spectral_sp_distribution <- rast(path)
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
Index_Alpha   = c('Shannon')
map_alpha_div(Input_Image_File = Input_Image_File,
              Output_Dir = Output_Dir,
              TypePCA = TypePCA,
              window_size = window_size,
              nbCPU = nbCPU,
              MaxRAM = MaxRAM,
              Index_Alpha = Index_Alpha,
              nbclusters = nbclusters)

# To have the name of the PCs selected to map the alpha diversity in the filename 
# run the function saved in the document map_alpha_div_PC_naming.R to have it in the 
# global environment

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
Shannon_map <- rast("~/data/biodivmapR_sent/RESULTS_aviris_extend/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_10_PC1278")

# Shannon_map <- rast("~/data/biodivmapR_sent/RESULTS/sent_crop_envi/SPCA/ALPHA/Shannon_10")
# Shannon_map <- rast("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_10")
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
  geom_spatraster(data = Shannon_map, na.rm = TRUE, aes(fill=Shannon_10_PC1278))+ #need to change the fill variable
  # theme_map()+
  theme_minimal()+
  scale_fill_gradientn(colours = viridis_colors, na.value = na_col) +
  labs(fill = substitute(paste("Shannon index", italic("H'"))))
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

sent_view <- rast("~/data/output/sent_crop_view.tif")
plotRGB(sent_view, r=3, g=2, b=1, scale=10000, stretch="lin", smooth=T)
sent_view_scale <-sent_view/7

boundary_strip_0708 <- st_read("~/scratch/reflectance_data/boundary_strip_0708.shp")
boundary_strip_0708_proj <- st_transform(boundary_strip_0708,32613)
area_interest <- st_read("areas_of_interest.gpkg")
area_interest_proj <- st_transform(area_interest,32613)
boundary_strip_0708_crop <- st_crop(boundary_strip_0708_proj, area_interest_proj)

s <- ggplot() +
  geom_spatraster_rgb(data = sent_view_scale, interpolate=T, r=3, g=2, b=1)+
  # geom_sf(data=boundary_strip_0708_crop, aes(alpha=0), show.legend = FALSE)+
  theme_map()
  # theme_minimal()
s

map_truecol <- s+
  ggspatial::annotation_scale(location="bl", pad_x=unit(0.5, "in"),
  pad_y = unit(0.4, "in"), style="ticks", line_col="white", text_col="white")+
  
  ggspatial::annotation_north_arrow(location="bl", which_north=T, pad_x=unit(0, "in"),
  pad_y = unit(0.4, "in"), height = unit(0.7, "cm"), width=unit(0.7, "cm"))
map_truecol

setwd("~/data/output")
save_plot(filename = "sentinel_full_view_aviris0708_boundary.png", map_truecol, base_asp = 2, bg = "white")
save_plot(filename = "sentinel_full_view.png", map_truecol, base_asp = 2, bg = "white")

cowplot::plot_grid(map_div, map_truecol, ncol = 2, align = "hv")


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

######################
# Beta diversity map
######################
print("MAP BETA DIVERSITY")
map_beta_div(Input_Image_File = Input_Image_File, 
             Output_Dir = Output_Dir, 
             TypePCA = TypePCA,
             window_size = window_size, 
             nbCPU = nbCPU, 
             MaxRAM = MaxRAM,
             nbclusters = nbclusters)

window_size <- 30
path <- paste(Output_Dir, "/", NameRaster, "/SPCA/BETA/BetaDiversity_BCdiss_PCO_", window_size, sep="")
Beta_map <- rast(path)
Beta_map_10 <- rast(path)
plotRGB(Beta_map, r=1, g=2, b=3, stretch="lin")

b <- ggplot() +
  geom_spatraster_rgb(data = Beta_map, r=1, g=2, b=3, interpolate = T)+
  # theme_map()+
  theme_minimal()
  # scale_fill_gradientn(colours = viridis_colors, na.value = na_col) +
  # labs(fill = substitute(paste("Shannon index", italic("H'"))))
b  


#######################################
# Comparison
#######################################
Shannon_PC12 <- rast("~/data/biodivmapR_sent/RESULTS_PC_selection/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_10_PC12")
Shannon_PC129 <- rast("~/data/biodivmapR_sent/RESULTS_PC_selection/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_10_PC129")
Shannon_PC12789 <- rast("~/data/biodivmapR_sent/RESULTS_PC_selection/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_10_PC12789")
par(mfrow=c(1,2))
plot(Shannon_PC12)
plot(Shannon_PC129)
plot(Shannon_PC12789)
par(mfrow=c(1,1))

viridis_colors <- viridis::plasma(20)
na_col <- grey(0.8, alpha=0.5)
map_shannon_PC12 <- ggplot() +
  geom_spatraster(data = Shannon_PC12, na.rm = TRUE, aes(fill=Shannon_10_PC12))+ #need to change the fill variable
  # theme_map()+
  theme_minimal()+
  scale_fill_gradientn(colours = viridis_colors, na.value = na_col) +
  labs(fill = substitute(paste("Shannon index", italic("H'"))))+
  plot_annotation(
    title = "Estimates of plant Shannon diversity index using PCs 1, 2",
    theme = theme(plot.title = element_text(hjust=0.4, size = 16, face = "bold")))+
  ggspatial::annotation_scale(location="bl", pad_x=unit(0.7, "in"),
                              pad_y = unit(0.6, "in"), style="ticks", line_col="black", text_col="black")+
  ggspatial::annotation_north_arrow(location="bl", which_north=T, pad_x=unit(0.7, "in"),
                                    pad_y = unit(0.8, "in"), height = unit(0.6, "cm"), width=unit(0.6, "cm"))
map_shannon_PC12  


map_shannon_PC129 <- ggplot() +
  geom_spatraster(data = Shannon_PC129, na.rm = TRUE, aes(fill=Shannon_10_PC129))+ #need to change the fill variable
  # theme_map()+
  theme_minimal()+
  scale_fill_gradientn(colours = viridis_colors, na.value = na_col) +
  labs(fill = substitute(paste("Shannon index", italic("H'"))))+
  plot_annotation(
    title = "Estimates of plant Shannon diversity index using PCs 1,2,9",
    theme = theme(plot.title = element_text(hjust=0.4, size = 16, face = "bold")))+
  ggspatial::annotation_scale(location="bl", pad_x=unit(0.7, "in"),
                              pad_y = unit(0.6, "in"), style="ticks", line_col="black", text_col="black")+
  ggspatial::annotation_north_arrow(location="bl", which_north=T, pad_x=unit(0.7, "in"),
                                    pad_y = unit(0.8, "in"), height = unit(0.6, "cm"), width=unit(0.6, "cm"))
map_shannon_PC129

  
cowplot::plot_grid(map_shannon_PC12, map_shannon_PC129, ncol = 2, align = "hv")
# install.packages("spatialEco")
library(spatialEco)
# install.packages("SpatialPack")
library(SpatialPack)
corr_diff_PC <- raster.modified.ttest(Shannon_PC12, Shannon_PC129)
plot(corr_diff_PC$corr, type="interval", breakby="cases")
plot(corr_diff_PC$p.value)
plot(corr_diff_PC$moran.x, type="interval", breakby="cases")
plot(corr_diff_PC$moran.y, type="interval", breakby="cases")
# writeRaster(corr_diff_PC, "~/data/biodivmapR_sent/RESULTS_PC_selection/sent_crop_envi_BIL/SPCA/ALPHA/correlation_PC12_PC129.envi")
corr_diff_PC <- rast("~/data/biodivmapR_sent/RESULTS_PC_selection/sent_crop_envi_BIL/SPCA/ALPHA/correlation_PC12_PC129.envi")
terra::plot(corr_diff_PC$corr, col=viridis_colors, type="continuous")
plot(corr_diff_PC$corr, type="interval", breakby="cases", main="correlation between Shannon index produced with PC 1,2 and PC 1,2,9")


map_shannon_PC12789 <- ggplot() +
  geom_spatraster(data = Shannon_PC12789, na.rm = TRUE, aes(fill=Shannon_10_PC12789))+ #need to change the fill variable
  # theme_map()+
  theme_minimal()+
  scale_fill_gradientn(colours = viridis_colors, na.value = na_col) +
  labs(fill = substitute(paste("Shannon index", italic("H'"))))+
  plot_annotation(
    title = "Estimates of plant Shannon diversity index using PCs 1,2,7,8,9",
    theme = theme(plot.title = element_text(hjust=0.4, size = 16, face = "bold")))+
  ggspatial::annotation_scale(location="bl", pad_x=unit(0.7, "in"),
                              pad_y = unit(0.6, "in"), style="ticks", line_col="black", text_col="black")+
  ggspatial::annotation_north_arrow(location="bl", which_north=T, pad_x=unit(0.7, "in"),
                                    pad_y = unit(0.8, "in"), height = unit(0.6, "cm"), width=unit(0.6, "cm"))
map_shannon_PC12789  
cowplot::plot_grid(map_shannon_PC12, map_shannon_PC129, map_shannon_PC12789, ncol = 2, align = "hv")
cowplot::plot_grid(map_shannon_PC12, map_shannon_PC12789, ncol = 2, align = "hv")



corr_diff_PC_12_12789 <- raster.modified.ttest(Shannon_PC12, Shannon_PC12789)
# writeRaster(corr_diff_PC_12_12789, "~/data/biodivmapR_sent/RESULTS_PC_selection/sent_crop_envi_BIL/SPCA/ALPHA/correlation_PC12_PC12789.envi")
corr_diff_PC_12_12789 <- rast("~/data/biodivmapR_sent/RESULTS_PC_selection/sent_crop_envi_BIL/SPCA/ALPHA/correlation_PC12_PC12789.envi")
plot(corr_diff_PC_12_12789$corr, type="interval", breakby="eqint", main="correlation between Shannon index produced with PC 1,2 and PC 1,2,7,8,9", col=viridis_colors)
plot(corr_diff_PC_12_12789$corr, type="interval", breaks=c(-1,-0.5,-0.1,0.1,0.5, 1), main="correlation between Shannon index produced with PC 1,2 and PC 1,2,7,8,9", col=viridis_colors[c(1,5,10,15,18)])
plot(corr_diff_PC_12_12789$corr, type="interval", breaks=c(-1,-0.5,-0.1,0.1,0.5, 1), main="correlation between Shannon index produced with PC 1,2 and PC 1,2,7,8,9", col=c(viridis_colors[c(4,8,12,16)], "#009200"), plg=list(title="correlation coefficient"))


############################
# Autocorrelation
############################
library(gstat)
library(sp)
library(raster)
Shannon_map <- rast("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_10")
Shannon_map_raster <- raster(Shannon_map)
Shannon_map_sp <- as(Shannon_map_raster, "SpatialPointsDataFrame")

vario_close <- variogram(Shannon_10~1, data=Shannon_map_sp, cutoff=5000, width=100)
plot(vario_close)
vario_far <- variogram(Shannon_10~1, data=Shannon_map_sp, cutoff=30000, width=1000)
plot(vario_far)

try <- vgm(0.1, "Mat", 2000, nugget=0.1)
plot(try, cutoff=15000)

vario_fit_close <- fit.variogram(vario_close, vgm(0.1, "Mat", 1000, nugget=0.1), fit.kappa = TRUE)
vario_fit_far <- fit.variogram(vario_far, vgm(0.05, "Mat", 4000, nugget=0.12), fit.kappa = TRUE)

plot(vario_fit_far, cutoff=5000)
plot(vario_far, vario_fit_far)


plot(vario_fit_close, cutoff=5000)
plot(vario_close, vario_fit_close)


##################
# try out distance
##################
Shannon <- rast(paste("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon", 10, "PC1278", sep="_"))
dist <- terra::distance(Shannon)
plot(dist)
Shannon
Shannon_for_dist <- ifel(is.na(Shannon), 999, NA)
dist_2 <- terra::distance(Shannon_for_dist)
plot(dist_2)
table(as.vector(is.na(dist_2)))
dist_2[dist_2 == 0] <- NA
table(as.vector(is.na(dist_2)))
table(as.vector(is.na(Shannon)))
plot(dist_2, colNA="red")
plot(Shannon, add=T, colNA="blue")
writeRaster(dist_2, filename = "~/data/biodivmapR_sent/dist_to_water.tif")


