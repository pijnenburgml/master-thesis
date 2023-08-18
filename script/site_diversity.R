setwd("~/data/field_work/")
library(sf)
library(terra)
library(raster)
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


########
# biodivmapR from plot
########
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
Kmeans_info <- get(load("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/SpectralSpecies/Kmeans_info.Rdata"))
PCA_Output <- get(load("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/PCA/PCA_Output.RData"))
# location of the spectral species raster needed for validation
Path_SpectralSpecies <- Kmeans_info$SpectralSpecies
nbclusters <- 20
# get diversity indicators corresponding to shapefiles (no partitioning of spectral dibversity based on field plots so far...)
Biodiv_Indicators <- diversity_from_plots_ML(Raster_SpectralSpecies = Path_SpectralSpecies, 
                                          Plots = Path_Vector,
                                          nbclusters = nbclusters)

Shannon_RS <- c(Biodiv_Indicators$Shannon)[[1]]

Spectral_sp_whole <- rast("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/SpectralSpecies/SpectralSpecies")
plot(Spectral_sp_whole$`Iter 1`)
site_1 <- vect("~/data/field_work/site_1.shp")
Spectral_sp_whole_crop <- crop(Spectral_sp_whole, site_1)
plot(Spectral_sp_whole_crop$`Iter 1`)
table(as.data.frame(Spectral_sp_whole_crop$`Iter 1`))


##############################
# Comparison with biodivmapR on the subset
##############################
site_1 <- vect("~/data/field_work/site_1.shp")
site_1_sent <- crop(sent_view, site_1)
plotRGB(site_1_sent, r=3, g=2, b=1, scale=1000, stretch="lin")
# writeRaster(site_1_sent, filename = "~/data/field_work/site_1_sent", filetype="ENVI", gdal=c("INTERLEAVE=BIL"), overwrite=T)
site_1_sent <- rast("~/data/field_work/site_1_sent")

mask <- rast("~/data/biodivmapR_sent/mask_sent2_final_NA")
mask_crop_site_1 <- crop(mask, site_1)
plot(mask_crop_site_1)
# writeRaster(mask_crop_site_1, filename = "~/data/field_work/mask_site_1", filetype="ENVI", gdal="INTERLEAVE=BSQ", overwrite=T, datatype="INT1U")
mask_site_1 <- rast("~/data/field_work/mask_site_1")

Datadir <- "~/data/field_work"
NameRaster <- "site_1_sent"
Input_Image_File <- file.path(Datadir,NameRaster)
NameMask <- "mask_site_1"
Input_Mask_File <- file.path(Datadir, NameMask)
Output_Dir <- '~/data/field_work/RESULTS'
# dir.create(path = Output_Dir,recursive = T,showWarnings = F)
Continuum_Removal <- F
TypePCA <- 'SPCA'
FilterPCA <- F
window_size <-30
# # computational parameters
nbCPU <- 2
MaxRAM <- 6
# number of clusters (spectral species)
nbclusters <- 10
print("PERFORM DIMENSIONALITY REDUCTION")
PCA_Output <- perform_PCA(Input_Image_File = Input_Image_File,
                          Input_Mask_File = Input_Mask_File,
                          Output_Dir = Output_Dir,
                          TypePCA = TypePCA,
                          FilterPCA = FilterPCA,
                          Continuum_Removal = Continuum_Removal,
                          nbCPU = nbCPU,
                          MaxRAM = MaxRAM)

Input_Mask_File <- PCA_Output$MaskPath

# way of visualizing PC variance explained
var_exp <- (PCA_Output$PCA_model$sdev^2/sum(PCA_Output$PCA_model$sdev^2))*100
barplot(var_exp, names.arg = colnames(PCA_Output$PCA_model$x))
# visual assesment of PCs
viridis_colors <- viridis::inferno(60)
PCs <- rast(file.path(Output_Dir, NameRaster, TypePCA, "PCA", "OutputPCA_10_PCs"))
plot(PCs, col=viridis_colors, legend=F)
SelectedPCs = c(1,2,4,7)

print("MAP SPECTRAL SPECIES")
future.seed=TRUE
Kmeans_info <- map_spectral_species(Input_Image_File = Input_Image_File,
                                    Input_Mask_File = PCA_Output$MaskPath,
                                    Output_Dir = Output_Dir,
                                    SpectralSpace_Output = PCA_Output,
                                    nbclusters = nbclusters,
                                    nbCPU = nbCPU, MaxRAM = MaxRAM,
                                    SelectedPCs = SelectedPCs)
Spectral_sp <- rast(file.path(Output_Dir, NameRaster, TypePCA, "SpectralSpecies","SpectralSpecies"))
table(as.data.frame(Spectral_sp$`Iter 1`))
plot(Spectral_sp$`Iter 1`)

Index_Alpha   = c('Shannon')
map_alpha_div(Input_Image_File = Input_Image_File,
              Output_Dir = Output_Dir,
              TypePCA = TypePCA,
              window_size = window_size,
              nbCPU = nbCPU,
              MaxRAM = MaxRAM,
              Index_Alpha = Index_Alpha,
              nbclusters = nbclusters)
# get error
# Error in alphaSSD[[i]] : subscript out of bounds 
Spectral_sp_df <- as.data.frame(Spectral_sp)
Spectral_sp_df_iter1 <- as.data.frame(Spectral_sp$`Iter 1`)
Spectral_species_vegan <- as.data.frame(matrix(ncol=10, nrow=ncol(Spectral_sp_df)))
colnames(Spectral_species_vegan) <- 1:10
rownames(Spectral_species_vegan) <- colnames(Spectral_sp_df)
table(as.data.frame(S))
Spectral_species_vegan[1,] <- as.numeric((table(Spectral_sp_df$`Iter 1`)))[-c(1)]
vegan::diversity(Spectral_species_vegan[1,])

#########################
# try with Aviris tile
#########################
Datadir <- "/scratch/mpijne/reflectance_data"
NameRaster <- "x_10_y_14"
# Define path for image file to be processed
Input_Image_File <- file.path(Datadir,NameRaster)
Output_Dir <- '~/scratch/biodivmapR_Aviris'
# dir.create(path = Output_Dir,recursive = T,showWarnings = F)
Continuum_Removal <- F
TypePCA <- 'SPCA'
FilterPCA <- F
window_size <-10
# # computational parameters
nbCPU <- 2
MaxRAM <- 6
# number of clusters (spectral species)
nbclusters <- 10
print("PERFORM DIMENSIONALITY REDUCTION")
PCA_Output <- perform_PCA(Input_Image_File = Input_Image_File,
                          Input_Mask_File = Input_Mask_File,
                          Output_Dir = Output_Dir,
                          TypePCA = TypePCA,
                          FilterPCA = FilterPCA,
                          Continuum_Removal = Continuum_Removal,
                          nbCPU = nbCPU,
                          MaxRAM = MaxRAM)

Input_Mask_File <- PCA_Output$MaskPath

# way of visualizing PC variance explained
var_exp <- (PCA_Output$PCA_model$sdev^2/sum(PCA_Output$PCA_model$sdev^2))*100
barplot(var_exp, names.arg = colnames(PCA_Output$PCA_model$x))
# visual assesment of PCs
viridis_colors <- viridis::inferno(60)
PCs <- rast(file.path(Output_Dir, NameRaster, TypePCA, "PCA", "OutputPCA_10_PCs"))
plot(PCs, col=viridis_colors, legend=F)
SelectedPCs = c(1,2,4,7)

print("MAP SPECTRAL SPECIES")
future.seed=TRUE
Kmeans_info <- map_spectral_species(Input_Image_File = Input_Image_File,
                                    Input_Mask_File = PCA_Output$MaskPath,
                                    Output_Dir = Output_Dir,
                                    SpectralSpace_Output = PCA_Output,
                                    nbclusters = nbclusters,
                                    nbCPU = nbCPU, MaxRAM = MaxRAM,
                                    SelectedPCs = SelectedPCs)









