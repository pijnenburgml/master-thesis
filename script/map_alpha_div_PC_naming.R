map_alpha_div_ML <- function(Input_Image_File = FALSE,
                          Input_Mask_File = FALSE,
                          Output_Dir = '',
                          window_size = 10,
                          TypePCA = "SPCA",
                          nbclusters = 50,
                          MinSun = 0.25,
                          pcelim = 0.02,
                          Index_Alpha = "Shannon",
                          FullRes = FALSE, LowRes = TRUE, MapSTD = TRUE,
                          nbCPU = 1, MaxRAM = 0.25, ClassifMap = FALSE, SelectedPCs=SelectedPCs) {
  
  # 1- get path for spectral species path, and possibly update Input_Image_File
  # and nbclusters if using classification map as input data
  SSDpathlist <- get_SSpath(Output_Dir, Input_Image_File, TypePCA, ClassifMap, nbclusters)
  Spectral_Species_Path <- SSDpathlist$Spectral_Species_Path
  PCs <- paste0("SpectralSpecies",paste0(SelectedPCs, collapse = ""))
  Spectral_Species_Path <-  file.path(Output_Dir_SS, PCs)
  SSD_Dir <- SSDpathlist$SSD_Dir
  Input_Image_File <- SSDpathlist$Input_Image_File
  nbclusters <- SSDpathlist$nbclusters
  
  # 2- COMPUTE ALPHA DIVERSITY
  ALPHA <- compute_alpha_metrics(Spectral_Species_Path = Spectral_Species_Path,
                                 SSD_Dir = SSD_Dir,
                                 window_size = window_size,
                                 Input_Mask_File = Input_Mask_File,
                                 nbclusters = nbclusters,
                                 MinSun = MinSun,
                                 pcelim = pcelim,
                                 nbCPU = nbCPU,
                                 MaxRAM = MaxRAM,
                                 Index_Alpha = Index_Alpha)
  
  # 3- SAVE ALPHA DIVERSITY MAPS
  print("Write alpha diversity maps")
  # which spectral indices will be computed
  Shannon <- Simpson <- Fisher <- FALSE
  if (length((grep("Shannon", Index_Alpha))) > 0) Shannon <- TRUE
  if (length((grep("Simpson", Index_Alpha))) > 0) Simpson <- TRUE
  if (length((grep("Fisher", Index_Alpha))) > 0) Fisher <- TRUE
  
  Output_Dir_Alpha <- define_output_subdir(Output_Dir, Input_Image_File, TypePCA, "ALPHA")
  HDR <- ALPHA$HDR
  if (Shannon == TRUE) {
    Index <- "Shannon"
    HDR$`band names` <- Index
    nb_PCs <- paste0(SelectedPCs, collapse = "")
    Alpha_Path <- file.path(Output_Dir_Alpha, paste(Index, "_", window_size, "_PC", nb_PCs, sep = ""))
    write_raster(Image = ALPHA$Shannon, HDR = HDR, ImagePath = Alpha_Path,
                 window_size = window_size, FullRes = FullRes, LowRes = LowRes,
                 SmoothImage = TRUE)
    if (MapSTD == TRUE) {
      Index <- "Shannon_SD"
      HDR$`band names` <- Index
      Alpha_Path <- file.path(Output_Dir_Alpha, paste(Index, "_", window_size, "_PC", nb_PCs, sep = ""))
      write_raster(ALPHA$Shannon_SD, HDR, Alpha_Path, window_size, FullRes = FullRes, LowRes = LowRes, SmoothImage = TRUE)
    }
  }
  
  if (Fisher == TRUE) {
    Index <- "Fisher"
    HDR$`band names` <- Index
    Alpha_Path <- file.path(Output_Dir_Alpha, paste(Index, "_", window_size, sep = ""))
    write_raster(ALPHA$Fisher, HDR, Alpha_Path, window_size, FullRes = FullRes, LowRes = LowRes, SmoothImage = TRUE)
    if (MapSTD == TRUE) {
      Index <- "Fisher_SD"
      HDR$`band names` <- Index
      Alpha_Path <- file.path(Output_Dir_Alpha, paste(Index, "_", window_size, sep = ""))
      write_raster(ALPHA$Fisher_SD, HDR, Alpha_Path, window_size, FullRes = FullRes, LowRes = LowRes, SmoothImage = TRUE)
    }
  }
  
  if (Simpson == TRUE) {
    Index <- "Simpson"
    HDR$`band names` <- Index
    Alpha_Path <- file.path(Output_Dir_Alpha, paste(Index, "_", window_size, sep = ""))
    write_raster(ALPHA$Simpson, HDR, Alpha_Path, window_size, FullRes = FullRes, LowRes = LowRes, SmoothImage = TRUE)
    if (MapSTD == TRUE) {
      Index <- "Simpson_SD"
      HDR$`band names` <- Index
      Alpha_Path <- file.path(Output_Dir_Alpha, paste(Index, "_", window_size, sep = ""))
      write_raster(ALPHA$Simpson_SD, HDR, Alpha_Path, window_size, FullRes = FullRes, LowRes = LowRes, SmoothImage = TRUE)
    }
  }
  return(invisible())
}

