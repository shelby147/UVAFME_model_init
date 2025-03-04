#
# #########################
# Purpose: Parameterize Roman
# Author: Adrianna C. Foster
# Date: August, 2022
# R version 4.2.0 (2022-04-22) 'Vigorous Calisthenics'
# #########################
# #########################
# Input format: csv, shp, tif
# Output format: csv
# #########################

path = ""
setwd(path)

library(dplyr)
library(data.table)
library(raster)
library(rgdal)
library(reshape2)
source('scripts/TreeLine_Funcs.R')

## parameterization files ------------------------------------------------------

climNA_fname <- 'raw/climate_dat_TVSF/climateNAinput_1961-1990M.csv'
cldmn_dir <- 'params_data/climate_dat/cld'
cldsd_dir <- 'params_data/climate_dat/cld_std'
windmn_dir <- 'params_data/climate_dat/wind/wc2.1_30s_wind'
lightning_dir <- 'params_data/climate_dat/lightning'
lightning_file <- 'AKLDN_10km'

## climate NA output
climNA_output <- fread(climNA_fname)

## site locations
# locs <- dplyr::select(climNA_output, id1, Latitude, Longitude) %>%
#   dplyr::rename('site' = 'id1', 'longitude' = 'Longitude',
#                 'latitude' = 'Latitude') %>%
#   distinct_all()

locs <- dplyr::select(climNA_output, ID1, Latitude, Longitude) %>%
  dplyr::rename('site' = 'ID1', 'longitude' = 'Longitude',
                'latitude' = 'Latitude') %>%
  distinct_all()


## parameterization ------------------------------------------------------------

## make climate data
clim_dat <- param_clim(climNA = climNA_output, cldmn_dir, cldsd_dir, 
                       windmn_dir)

## filter out sites with missing climate data
locs <- filter(locs, !(site %in% clim_dat$bad_sites))

## make lightning strike data
lightning_dat <- param_lightning(lightning_dir, lightning_file, locs)

## write to files --------------------------------------------------------------

outdir <- 'temp'

write.csv(clim_dat$clim_mn, paste0(outdir, '/UVAFME2018_climate.csv'),
          row.names = F, quote = F, na = '')
write.csv(clim_dat$clim_sd, paste0(outdir, '/UVAFME2018_climate_stddev.csv'),
          row.names = F, quote = F, na = '')
write.csv(clim_dat$extra_mn, paste0(outdir, '/UVAFME2018_climate_ex.csv'),
          row.names = F, quote = F, na = '')
write.csv(clim_dat$extra_sd, paste0(outdir, '/UVAFME2018_climate_ex_stddev.csv'),
          row.names = F, quote = F, na = '')
write.csv(lightning_dat, paste0(outdir, '/UVAFME2018_lightning.csv'),
          row.names = F, quote = F, na = '')
