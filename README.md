# Mapping of plant diversity in the Arctic tundra using a spectral diversity measure based on the spectral species concept

## Content
This repository contains all the code to produce the study presented in the master thesis "Mapping of spectral species diversity in the Arctic tundra and its relationship with topographic complexity – a use case of the spectral species concept in biodiversity research and conservation", done in 2023 by Marie-Lou Pijnenburg (marielou.pijnenburg2000@gmail.com)

## Data
* The Sentinel-2 data can be download from the open hub from copernicus schihub (https://scihub.copernicus.eu/). The two tile used are S2B_MSIL2A_20190727T185929_N0213_R013_T13WES_20190727T214238.SAFE and S2B_MSIL1C_20190727T185929_N0208_R013_T13WDS_20190727T210028.SAFE

* The data from the ABoVE campaign can be found and download from the earth data website under the NASA project "ABoVE: Hyperspectral Imagery AVIRIS-NG, Alaskan and Canadian Arctic, 2017-2019 V2" (Miller et al., 2022). The flight strip use has the following ID: ang20190802t220708.

* The elevation data can be found and download from the Arctic Digital Elevation Model (Porter et al., 2023) website https://www.pgc.umn.edu/data/arcticdem. We use mosaic tile at 2 meters resolution, from the version 4.1. The tile IDs are “29_21_1_” and “29_21_2_1”. 

* The files area_of_interest.gpkg, site_boundaries.gpkg and ang20190802t220708_outline_KML.kml can be found in this directory

## Data analysis

The data analysis and the production of the plots can be done with the code found under the script directory of this repository, following this order:

* The R scripts with name "map_spectral_sp_PC_naming", "map_alpha_div_PC_naming" and "hatched_function" are functions written by different authors, to which I included small changes to adapt them to this study. 
* The script Sentinel2019_preparation prepare the Sentinel-2 data and the corresponding mask
* The script BiodivmapR_Sentinel run the workflow from the biodivmapR package on the Sentinel-2 data
* The script site_diversity get the spectral diversity from the fieldwork sites
* The script rf_cloud_detection builds the RandomForest classifier for detection of cloudy pixels in the AVIRIS-NG data
* The script Aviris_data_preparation prepare the data from the ABoVE campaign and produce the corresponding mask
* The script biodivmapR_aviris run the biodivmapR workflow on the data from the ABoVE campaign 
* The script ArcDEM_preparation prepare the data from the Arctic Digital Elevation Model
* The script scaling_analysis contain the code to get the elevation data as well as the spectral species diversity data derived from Sentinel-2 data at different scales and the code to run the corresponding statistical models. 
* The script aviris_modelling contains the code to get the elevation data as well as the specral sepcies diversity data derived from the AVIRIS-NG data at different scales and the cade to run the corresponding statistical models. 
* The script Sentinel2_mapping contains the code to produce the Figure 1 of this project
* The script semi_variogram_plot contains the code to generates the plot of semi-variogram found in the Supplementary material
* The script fieldwork_table.Rmd produce the Supplementary Table S1 and require the a dataset "presence_absence_cleaned.csv" that is currently not publicly available but can be obtained via Jakob J. Assmann (jakob.assmann@uzh.ch)
* The script table_formation.Rmd produce the supplementary tables S2-S5. 

## License
This work is under the MIT license, which can be found under the text file called "license"

## Reference
Miller, C. E., Green, R. O., Thompson, D. R., Thorpe, A. K., Eastwood, M., Mccubbin, I. B., Olson-Duvall, W., Bernas, M., Sarture, C. M., Nolte, S., Rios, L. M., Hernandez, M. A., Bue, B. D., & Lundeen, S. R. (2022). ABoVE: Hyperspectral imagery AVIRIS-NG, Alaskan and Canadian Arctic, 2017-2019 V2. ORNL Distributed Active Archive Center. https://doi.org/10.3334/ORNLDAAC/2009

Porter, C., Howat, I., Noh, M.-J., Husby, E., Khuvis, S., Danish, E., Tomko, K., Gardiner, J., Negrete, A., Yadav, B., Klassen, J., Kelleher, C., Cloutier, M., Bakker, J., Enos, J., Arnold, G., Bauer, G., & Morin, P. (2023). ArcticDEM - Mosaics, Version 4.1. Harvard Dataverse. https://doi.org/10.7910/DVN/3VDC4W





