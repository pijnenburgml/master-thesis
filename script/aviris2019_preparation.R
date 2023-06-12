setwd("~/scratch/")
library(sf)
library(tidyverse)
library(dplyr)
library(terra)
library(caTools)
######
# 
######

line47 <- read.ENVI("reflectance_data/ang20190801t160747_rfl")
