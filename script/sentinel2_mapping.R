setwd("~/data")
library(sf)
library(tidyverse)
library(dplyr)
library(terra)
# library(biodivMapR)

############################
# area of interest
############################
area_interest <- st_read("areas_of_interest.gpkg")
area_interest_proj <- st_transform(area_interest,32613)

###########
# maps overview of the world
###########
library(maptools)
library(sp)
library(rgdal)
require(maps)
library(ggplot2)
library(cowplot)

lim = data.frame(ylim=c(25, 75), xlim=c(-150, -70))
bbox_us_ca <- st_bbox(c(
  xmin=lim$xlim[1],
  xmax=lim$xlim[2],
  ymin=lim$ylim[1],
  ymax=lim$ylim[2]))

area_interest_point <- data.frame(x=(st_bbox(area_interest)[1]+st_bbox(area_interest)[3])/2, y=(st_bbox(area_interest)[2]+st_bbox(area_interest)[4])/2)

map_us_ca <- ggplot()+
  borders("world", fill="grey90",colour="grey")+
  # coord_fixed(ylim=lim$ylim, xlim=lim$xlim)+
  coord_map(projection = "mercator", xlim=lim$xlim, ylim=lim$ylim, orientation = c(90, 0, 0))+
  geom_point(data=area_interest_point, aes(x,y, shape=15))+
  scale_shape_identity()+
  theme_map()
map_us_ca
save_plot(map_us_ca, filename = "~/data/output/final_plot/overview_us_ca_map.png", base_height = 3, bg="white")


# Plotting of RGB view

sent_view <- rast("~/data/output/sent_crop_view.tif")
plotRGB(sent_view, r=3, g=2, b=1, scale=10000, stretch="lin", smooth=T)
sent_view_scale <-sent_view/6
boundary_strip_0708 <- st_read("~/scratch/reflectance_data/boundary_strip_0708.shp")
boundary_strip_0708_proj <- st_transform(boundary_strip_0708,32613)
area_interest <- st_read("areas_of_interest.gpkg")
area_interest_proj <- st_transform(area_interest,32613)
boundary_strip_0708_crop <- st_crop(boundary_strip_0708_proj, area_interest_proj)
site_boundaries <-st_read("~/data/field_work/site_boundaries.gpkg")
site_boundaries$site <- c("Site 1", "Site 2", "Site 3")
bounding_box <- st_as_sf(as.polygons(ext(sent_view_scale)))
st_crs(bounding_box) <- crs(site_boundaries)

label_fill <- grey(0.9, alpha=0.5)
s <- ggplot() +
  geom_spatraster_rgb(data = sent_view_scale, interpolate=T, r=3, g=2, b=1)+
  geom_sf(data=boundary_strip_0708_crop, fill="white", alpha=0.5, colour="red", show.legend = FALSE)+
  geom_sf(data=site_boundaries, colour="red", fill="red", cex=5, show.legend = F)+
  geom_sf(data=bounding_box, colour="red", fill="white", alpha=0, linewidth=0.65)+
  geom_sf_label(data=site_boundaries, mapping=aes(label=site), fill="white", alpha=0.8, nudge_x = c(0), nudge_y = c(-2000),
                inherit.aes = F, label.size=0, show.legend = F)+
  geom_sf_label(boundary_strip_0708_crop, mapping=aes(label=paste("ABoVE campaign", "flight strip",sep="\n")), fill="white", alpha=0.8, nudge_x = c(-14000),
                nudge_y = c(12500), inherit.aes = F, label.size=0, show.legend = F)+
  geom_sf_label(data=bounding_box, mapping = aes(label="Sentinel-2 data"), fill="white", alpha=0.8, 
                nudge_y = c(12000), nudge_x = c(15000))+
  theme_map()
s

map_truecol <- s+
  ggspatial::annotation_scale(location="bl", pad_x=unit(0.5, "in"),
                              pad_y = unit(0.4, "in"), style="ticks", line_col="white", text_col="white")+
  
  ggspatial::annotation_north_arrow(location="bl", which_north=T, pad_x=unit(0, "in"),
                                    pad_y = unit(0.4, "in"), height = unit(0.7, "cm"), width=unit(0.7, "cm"))
map_truecol
setwd("~/data/output/final_plot")
# save_plot(filename = "sentinel_full_view_aviris0708_boundary_field_site.png", map_truecol, base_height = 4, bg = "white")
save_plot(filename = "sentinel_full_view_aviris0708_boundary_field_site.svg", map_truecol, base_height = 4, bg = "white")


