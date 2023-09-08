map_spectral_species_ML <- function(Input_Image_File, Output_Dir,
                                 SpectralSpace_Output,
                                 Input_Mask_File = FALSE,
                                 nbclusters = 50,
                                 nbCPU = 1, MaxRAM = 0.25,
                                 Kmeans_Only = FALSE, SelectedPCs = FALSE,
                                 SpectralFilter = NULL) {
  
  # check if input mask file has expected format
  if (!Input_Mask_File==FALSE){
    # driverMask_class <- new("GDALReadOnlyDataset", Input_Mask_File)
    driverMask_class <- rgdal::GDAL.open(filename = Input_Mask_File,
                                  read.only = TRUE, silent=FALSE,
                                  allowedDrivers = NULL, options=NULL)
    driverMask <- rgdal::getDriverLongName(rgdal::getDriver(driverMask_class))
    rgdal::GDAL.close(driverMask_class)
    if (driverMask == 'ENVI .hdr Labelled'){
      HDR <- read_ENVI_header(get_HDR_name(Input_Mask_File))
      if (!HDR$`data type`==1){
        Input_Mask_File <- check_update_mask_format(Input_Mask_File, Input_Image_File)
      }
    } else {
      Input_Mask_File <- check_update_mask_format(Input_Mask_File, Input_Image_File)
    }
  } else {
    message('Input_Mask_File not provided in function map_spectral_species.')
    message('Assuming all pixels are valid')
    message('A blank mask will be created for the need of next processing steps')
    
    HDR <- read_ENVI_header(get_HDR_name(Input_Image_File))
    Mask <- t(matrix(as.integer(1), nrow = HDR$lines, ncol = HDR$samples))
    MaskPath_Update <- paste(file_path_sans_ext(Input_Image_File),'_BlankMask',sep = '')
    Input_Mask_File <- update_shademask(MaskPath = FALSE,
                                        HDR = HDR,
                                        Mask = Mask,
                                        MaskPath_Update = MaskPath_Update)
  }
  
  Kmeans_info <- NULL
  # if no prior diversity map has been produced --> need PCA file
  if (is.null(SpectralSpace_Output$PCA_Files) | is.null(SpectralSpace_Output$TypePCA)){
    message('Please define input variable SpectralSpace_Output as a list including')
    message('PCA_Files: corresponds to the raster data to be processed (not necessarily resulting from PCA)')
    message('TypePCA: defines main directory where outputs will be written')
    message('This variable is automatically produced as an output of function perform_PCA()')
    message('However, you can set it manually, for example if you want to use spectral indices')
    message('as input raster data instead of PCA file produced from reflectance data')
    stop()
  }
  if (!file.exists(SpectralSpace_Output$PCA_Files)) {
    error_no_PCA_file(SpectralSpace_Output$PCA_Files)
    stop()
  }
  
  # define directories
  Output_Dir_SS <- define_output_subdir(Output_Dir, Input_Image_File, SpectralSpace_Output$TypePCA, "SpectralSpecies")
  Output_Dir_PCA <- define_output_subdir(Output_Dir, Input_Image_File, SpectralSpace_Output$TypePCA, "PCA")
  PCs <- paste0("SpectralSpecies",paste0(SelectedPCs, collapse = ""))
  Spectral_Species_Path <-  file.path(Output_Dir_SS, PCs)
  
  # 1- Select components used to perform clustering
  if (typeof(SelectedPCs) == 'logical'){
    if (SelectedPCs == FALSE){
      PC_Select_Path <- file.path(Output_Dir_PCA, "Selected_Components.txt")
    } else {
      message('Error when defining SelectedPCs :')
      message('either set SelectedPCs = FALSE')
      message('or provide a vectorincluding the rank of the variables to be selected from SpectralSpace_Output$PCA_Files')
      stop()
    }
  } else {
    PC_Select_Path = 'NoFile'
  }
  if (file.exists(PC_Select_Path)) {
    PC_Select <- utils::read.table(PC_Select_Path)[[1]]
  } else if (is.numeric(SelectedPCs)){
    PC_Select <- SelectedPCs
  } else {
    error_PC_sel(Output_Dir_PCA)
    stop()
  }
  message("Selected components:")
  print(PC_Select)
  
  # 2- sample data from PCA image
  ImNames <- list(Input_Image = Input_Image_File,
                  Mask_list = Input_Mask_File)
  if (is.null(SpectralSpace_Output$nb_partitions)){
    nb_partitions <- 20
  } else {
    nb_partitions <- SpectralSpace_Output$nb_partitions
  }
  Pix_Per_Partition <- define_pixels_per_iter(ImNames, nb_partitions = nb_partitions)
  
  ImPathHDR <- get_HDR_name(SpectralSpace_Output$PCA_Files)
  HDR <- read_ENVI_header(ImPathHDR)
  Subset <- get_random_subset_from_image(ImPath = SpectralSpace_Output$PCA_Files,
                                         MaskPath = Input_Mask_File,
                                         nb_partitions = nb_partitions,
                                         Pix_Per_Partition = Pix_Per_Partition,
                                         kernel = NULL,MaxRAM = MaxRAM)
  SubsetInit <- Subset
  dataPCA <- Subset$DataSubset[, PC_Select]
  if (length(PC_Select) == 1) {
    dataPCA <- matrix(dataPCA, ncol = 1)
  }
  
  # 3- PERFORM KMEANS FOR EACH ITERATION & DEFINE SPECTRAL SPECIES
  print("perform k-means clustering for each subset and define centroids")
  Kmeans_info <- init_kmeans(dataPCA = dataPCA,
                             nb_partitions = nb_partitions,
                             nbclusters = nbclusters,
                             nbCPU = nbCPU)
  Kmeans_info$SpectralSpecies <- Spectral_Species_Path
  
  if (Kmeans_info$Error==FALSE){
    if (Kmeans_Only==FALSE){
      ##    3- ASSIGN SPECTRAL SPECIES TO EACH PIXEL
      apply_kmeans(PCA_Path = SpectralSpace_Output$PCA_Files,
                   PC_Select = PC_Select,
                   Input_Mask_File = Input_Mask_File,
                   Kmeans_info = Kmeans_info,
                   Spectral_Species_Path = Spectral_Species_Path,
                   nbCPU = nbCPU, MaxRAM = MaxRAM)
    } else {
      print("'Kmeans_Only' was set to TRUE: kmeans was not applied on the full image")
      print("Please set 'Kmeans_Only' to FALSE if you want to produce spectral species map")
    }
    # save kmeans info into binary variable
    Kmeans_Path <- file.path(Output_Dir_PCA, "Kmeans_Info.RData")
    save(Kmeans_info, file = Kmeans_Path)
  } else {
    ##    produce error report
    # create directory where error should be stored
    Output_Dir_Error <- define_output_subdir(Output_Dir, Input_Image_File, SpectralSpace_Output$TypePCA, "ErrorReport")
    # identify which samples cause problems
    LocError <- unique(c(which(!is.finite(Kmeans_info$MinVal)),which(!is.finite(Kmeans_info$MaxVal))))
    ValError <- which(!is.finite(dataPCA[,LocError[1]]))
    # Get the original data corresponding to the first sample
    DataError <- SubsetInit$DataSubset[ValError,]
    DataErrorCR <- Subset$DataSubset[ValError,]
    CoordinatesError <- SubsetInit$coordPix[ValError,]
    # save these in a RData file
    FileError <- file.path(Output_Dir_Error,'ErrorPixels.RData')
    ErrorReport <- list('CoordinatesError' = CoordinatesError,'DataError' = DataError,
                        'DataError_afterCR' = DataErrorCR, 'SpectralFilter'=SpectralFilter)
    save(ErrorReport, file = FileError)
    message("")
    message("*********************************************************")
    message("       An error report directory has been produced.      ")
    message("Please check information about data causing errors here: ")
    print(FileError)
    message("               The process will now stop                 ")
    message("*********************************************************")
    message("")
    stop()
  }
  return(Kmeans_info)
}
