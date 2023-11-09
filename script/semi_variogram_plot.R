############################
# Plot autocorrelation
############################
library(gstat)
library(sp)
library(raster)
Shannon_map <- rast("~/data/biodivmapR_sent/RESULTS_cluster_20/sent_crop_envi_BIL/SPCA/ALPHA/Shannon_10_PC1278")
Shannon_map_raster <- raster(Shannon_map)
Shannon_map_sp <- as(Shannon_map_raster, "SpatialPointsDataFrame")
vario_close <- variogram(Shannon_10_PC1278~1, data=Shannon_map_sp, cutoff=5000, width=100)
plot(vario_close)
try <- vgm(0.1, "Mat", 500, nugget=0.1)
plot(try, cutoff=5000)
vario_fit_close <- fit.variogram(vario_close, try, fit.kappa = TRUE)
plot(vario_fit_close, cutoff=5000)
plot(vario_close, vario_fit_close)
vario_fit_close
range_sh <- round(vario_fit_close$range[2], digits = 2)
preds = variogramLine(vario_fit_close, maxdist = 5000)
semivar_shannon <- ggplot()+
  geom_point(data=vario_close, aes(x=dist, y=gamma), pch=21, cex=2,  col="dodgerblue1")+
  geom_line(data=preds, aes(x=dist, y=gamma), col="dodgerblue1")+
  labs(
    y=expression(atop("semi-variance of Sentinel-2-derived", 
                      "Shannon index")), 
    x=("distance [m]"))+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), text=element_text(size=14),
        axis.title.y = element_text(vjust=3), 
        plot.margin=margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"))+
  coord_cartesian(xlim =c(250, 5000), ylim = c(0, 0.15))+
  annotate("text", x=3500, y=0.02, label = paste("Range:", range_sh, sep=" "), size=6)
semivar_shannon

# save_plot(plot, filename = "~/data/output/final_plot/variogram_close_shannon_res100_width100_cutoff5000.png", base_height = 4)

# sd (elev)
sd_ele_100 <- rast("~/data/ArcDEM/ArcDEM_masked_10_res.tif")
sd_ele_100_raster <- raster(sd_ele_100)
sd_ele_100_sp <- as(sd_ele_100_raster, "SpatialPointsDataFrame")
vario_log_ele_close <- variogram(log(X29_21_1_1_2m_v4.1_dem)~1, data=sd_ele_100_sp, cutoff=5000, width=100)
plot(vario_log_ele_close)
try <- vgm(0.1, "Exp", 1000, nugget=0.12)
plot(try, cutoff=5000)
vario_fit_close <- fit.variogram(vario_log_ele_close, try, fit.kappa = TRUE)
plot(vario_log_ele_close, vario_fit_close, xlab="distance [m]")
vario_fit_close
plot(vario_log_ele_close)
plot(vario_fit_close, cutoff=5000)
plot(vario_log_ele_close, vario_fit_close)
range_ele <- round(vario_fit_close$range[2], digits = 2)
preds = variogramLine(vario_fit_close, maxdist = 5000)
semivar_ele <- ggplot()+
  geom_point(data=vario_log_ele_close, aes(x=dist, y=gamma), pch=21, cex=2,  col="dodgerblue1")+
  geom_line(data=preds, aes(x=dist, y=gamma), col="dodgerblue1")+
  labs(y=expression(paste("semi-variance of ", sigma, " (elevation)")), x=("distance [m]"))+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), text=element_text(size=14),
        plot.margin=margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"))+
  coord_cartesian(xlim =c(250, 5000), ylim = c(0, 0.7))+
  annotate("text", x=3500, y=0.1, label = paste("Range:", range_ele, sep=" "), size=6)
semivar_ele

# sd(slope)
sd_slope_100 <- rast("~/data/ArcDEM/sd_slope_masked_10_res.tif")
sd_slope_100_raster <- raster(sd_slope_100)
sd_slope_100_sp <- as(sd_slope_100_raster, "SpatialPointsDataFrame")
vario_slope_close <- variogram(log(slope)~1, data=sd_slope_100_sp, cutoff=5000, width=100)
plot(vario_slope_close)
try <- vgm(0.1, "Exp", 1000, nugget=0.12)
plot(try, cutoff=5000)
vario_fit_slope_close <- fit.variogram(vario_slope_close, try, fit.kappa = TRUE)
plot(vario_slope_close, vario_fit_slope_close, xlab="distance [m]")
vario_fit_slope_close
range_sp <- round(vario_fit_slope_close$range[2], digits = 2)

preds = variogramLine(vario_fit_slope_close, maxdist = 5000)
semivar_slope <- ggplot()+
  geom_point(data=vario_slope_close, aes(x=dist, y=gamma), pch=21, cex=2,  col="dodgerblue1")+
  geom_line(data=preds, aes(x=dist, y=gamma), col="dodgerblue1")+
  labs(y=expression(paste("semi-variance of ", sigma, " (slope)")), x=("distance [m]"))+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), text=element_text(size=14), 
        plot.margin=margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"))+
  coord_cartesian(xlim =c(250, 5000), ylim = c(0, 0.4))+
  annotate("text", x=3500, y=0.05, label = paste("Range:", range_sp, sep=" "), size=6)
semivar_slope

vario_together <- plot_grid(semivar_shannon, semivar_slope, semivar_ele, align = "v", ncol=1, axis = "r",labels="AUTO", label_size = 25)
save_plot(vario_together, filename = "~/data/output/final_plot/variogram_ele_slope_shannon.png",base_height = 15, base_width = 6, bg="white")

