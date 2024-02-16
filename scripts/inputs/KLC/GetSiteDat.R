#
# #########################
# Purpose: Grab data for creating UVAFME2018_site.csv file
# Author: Adrianna C. Foster
# Date: July, 2020
# R version 4.0.2 (2020-06-22) 'Taking Off Again'
# #########################
# #########################
# Input format: csv, tiff, shp
# Output format: csv
# #########################

library(dplyr)
library(reshape2)
library(rgdal)
library(raster)

## inventory locations with siteID, name, lat/lon, stand_age
# locs <- read.csv('intermediate_data/sitedat_input.csv') %>%
#   dplyr::select(site, name, latitude, longitude)

locs <- read.csv('scripts/inputs/KLC/intermediate_data/Landsat_init_sitedat_input.csv') %>%
  dplyr::select(site, name, latitude, longitude)

#1deg_grid for grabbing correct ABoVE DEM tiles
latlon_grid <- readOGR(dsn = 'UVAFME_Parameterization/GIS_Data',
                       layer = 'above_latlon_1deg')
latlon_grid$left <- abs(round(latlon_grid$left, 0))
latlon_grid$bottom <- abs(round(latlon_grid$bottom, 0))
latlon_grid$id <- paste0('N', latlon_grid$bottom, 'W',
                         stringr::str_pad(latlon_grid$left, width = 3,
                                          side = 'left', pad = '0'))

capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}



get_siteDat <- function(locs) {
  ## function to parameterize sites from input topography rasters
  ## Inputs:
  ##  locs - site locations, should have site, name, latitude, longitude, stand_age

  ## convert to spdf
  wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  locs_spdf <- SpatialPointsDataFrame(coords = locs[,c('longitude','latitude')],
                                      data = locs,
                                      proj4string = wgs)

  ## get 1-degree grid location for each site
  ex_latlon <- over(spTransform(locs_spdf, crs(latlon_grid)), latlon_grid)$id
  locs_spdf@data <- cbind(locs_spdf@data, ex_latlon)

  ## what are unique values for lat-lon grid locations
  latlons <- unique(ex_latlon)

  ## grab data from ABoVE DEM
  out <- data.frame()
  for (i in 1:length(latlons)){

    locs_spdf@data <- dplyr::select(locs_spdf@data, site:ex_latlon)
    sub_locs <- subset(locs_spdf, locs_spdf$ex_latlon %in% latlons[i])

    ## read in/calculate topography rasters
    dem <- raster(paste0('/projects/above_gedi/users/pburns/ABoVE',
                         '/DEM_mosaics/v2_1/all_ABoVE/2_DEMs/',
                         'ABoVE_composDEM_0p33AS_', latlons[i], '.tif'))

    slope <- raster(paste0('/projects/geode_data/share/elev/abovedem',
                           '/v2p1/slope_deg/ABoVE_composDEM_0p33AS_',
                           latlons[i], '_slope_deg.tif'))

    aspect <- raster(paste0('/projects/geode_data/share/elev/abovedem',
                           '/v2p1/aspect_deg/ABoVE_composDEM_0p33AS_',
                           latlons[i], '_aspect_deg.tif'))


    ## extract and add to locs spdf
    elevation <- extract(dem, spTransform(sub_locs, crs(dem)))
    slpe <- extract(slope, spTransform(sub_locs, crs(slope)))
    aspct <- extract(aspect, spTransform(sub_locs, crs(aspect)))
    sub_locs@data <- cbind(sub_locs@data, elevation, slpe, aspct)
    sub  <- as.data.frame(subset(sub_locs, !is.na(elevation) & !is.na(slpe) &
                                   !is.na(aspct)))
    out <- rbind(out, sub)
  }
  out <- dplyr::select(out, site:aspct)


  ## Create new spdf
  locs_spdf <- SpatialPointsDataFrame(coords = out[,c('longitude','latitude')],
                                      data = dplyr::select(out, site:aspct),
                                      proj4string = wgs)

  ## Read in soil maps
  soil_dir <- paste0('/projects/geode_data/share/geol_soil/',
                     'OpenLandMap/0_orig/')
  obd <- raster(paste0(soil_dir, 'sol_bulkdens.fineearth_usda.',
                       '4a1h_m_250m_b0..0cm_1950..2017_v0.2.tif'))
  abd <- raster(paste0(soil_dir, 'sol_bulkdens.fineearth_usda.',
                       '4a1h_m_250m_b30..30cm_1950..2017_v0.2.tif'))
  afc <- raster(paste0(soil_dir, 'sol_watercontent.33kPa_usda.4b1c_m_250m_',
                       'b30..30cm_1950..2017_v0.1.tif'))
  ofc <- raster(paste0(soil_dir, 'sol_watercontent.33kPa_usda.4b1c_md_250m_b',
                       '0..0cm_1950..2017_v0.1.tif'))
  apwp <- raster(paste0(soil_dir, 'sol_watercontent.1500kPa_usda.3c2a1a_m_250m_',
                        'b30..30cm_1950..2017_v0.1.tif'))
  opwp <- raster(paste0(soil_dir, 'sol_watercontent.1500kPa_usda.3c2a1a_m_250m_b',
                        '0..0cm_1950..2017_v0.1.tif'))
  sand <- raster(paste0(soil_dir, 'sol_sand.wfraction_usda.3a1a1a_m_250m_',
                        'b30..30cm_1950..2017_v0.2.tif'))
  twi <- raster(paste0(soil_dir, 'dtm_twi_merit.dem_m_500m_s0..0cm_2017_v1.0.tif'))

  ## Ecoregions
  eco <- readOGR(dsn = 'polygons/', layer = 'ABOVE_Ecoregions_L3')
  ## Extract all values
  ex_eco <- over(spTransform(locs_spdf, crs(eco)), eco)
  ex_obd <- extract(obd, spTransform(locs_spdf, crs(obd)))
  ex_abd <- extract(abd, spTransform(locs_spdf, crs(abd)))
  ex_afc <- extract(afc, spTransform(locs_spdf, crs(afc)))
  ex_ofc <- extract(ofc, spTransform(locs_spdf, crs(ofc)))
  ex_apwp <- extract(apwp, spTransform(locs_spdf, crs(apwp)))
  ex_opwp <- extract(opwp, spTransform(locs_spdf, crs(opwp)))
  ex_sand <- extract(sand, spTransform(locs_spdf, crs(sand)))
  ex_twi <- extract(twi, spTransform(locs_spdf, crs(twi)))

  ## Bind to data
  locs_spdf@data <- cbind(locs_spdf@data, ex_obd, ex_abd, ex_afc, ex_ofc,
                          ex_apwp, ex_opwp, ex_sand, ex_twi, ex_eco$NA_L2NAME)

  ## convert back to df, rename some columns
  siteDat <- as.data.frame(locs_spdf) %>%
    rename('slope' = 'slpe', 'aspect' = 'aspct', 'region' = 'ex_eco.NA_L2NAME')

  ## get texture
  siteDat$itxt <- NA
  siteDat$itxt[siteDat$ex_sand > 40] <- 1
  siteDat$itxt[siteDat$ex_sand <= 40] <- 2

  ## Convert units
  siteDat$a_fc <- siteDat$ex_afc/100.0
  siteDat$a_pwp <- siteDat$ex_apwp/100.0
  siteDat$o_fc <- siteDat$ex_ofc/100.0
  siteDat$o_pwp <- siteDat$ex_opwp/100.0
  siteDat$a_bd <- siteDat$ex_abd*10
  siteDat$o_bd <- siteDat$ex_obd*10

  ## We need saturation capacity
  siteDat$a_sat <- siteDat$a_pwp + siteDat$a_fc
  siteDat$o_sat <- siteDat$o_pwp + siteDat$o_fc

  siteDat$region <- capwords(tolower(siteDat$region))
  siteDat$region <- gsub(' ', '_', siteDat$region)

  ## Default values - can change later
  siteDat$hum_init <- 55.3
  siteDat$wprob <- 0.0

  ## bad data for a_sat, etc.
  bad_sites <- filter(siteDat, a_sat <= 0.0 | a_fc <= 0.0 |
                        a_pwp <= 0.0 | o_sat <= 0.0 | o_fc <= 0.0 |
                        o_pwp <= 0.0 | o_sat <= o_fc | o_fc <= o_pwp |
                        a_sat <= a_fc | a_fc <= a_pwp)
  bad_sites = bad_sites$site

  ## Use defaults
  siteDat$a_sat[siteDat$site %in% bad_sites] <- 0.44
  siteDat$a_fc[siteDat$site %in% bad_sites] <- 0.29
  siteDat$a_pwp[siteDat$site %in% bad_sites] <- 0.06

  siteDat$o_sat[siteDat$site %in% bad_sites] <- 0.38
  siteDat$o_fc[siteDat$site %in% bad_sites] <- 0.3
  siteDat$o_pwp[siteDat$site %in% bad_sites] <- 0.038


  ## Normalize TWI
  siteDat$norm_twi <- siteDat$ex_twi/132

  ## Get landscape position type
  siteDat$type <- NA
  siteDat$type[siteDat$slope <= 2.0 & siteDat$norm_twi > 0.75] <- 'Flat with Accumulation'
  siteDat$type[siteDat$slope <= 2.0 & siteDat$norm_twi <= 0.75] <- 'Flat without Accumulation'
  siteDat$type[siteDat$slope > 2.0 & siteDat$norm_twi > 0.75] <- 'Sloped with Accumulation'
  siteDat$type[siteDat$slope > 2.0 & siteDat$norm_twi <= 0.75] <- 'Sloped without Accumulation'

  ## Reset some values based on landscape position
  siteDat$a_sat[siteDat$type == 'Flat with Accumulation'] <- 0.53
  siteDat$a_fc[siteDat$type == 'Flat with Accumulation'] <- 0.38
  siteDat$a_pwp[siteDat$type == 'Flat with Accumulation'] <- 0.06
  siteDat$o_sat[siteDat$type == 'Flat with Accumulation'] <- 0.53
  siteDat$o_fc[siteDat$type == 'Flat with Accumulation'] <- 0.38
  siteDat$o_pwp[siteDat$type == 'Flat with Accumulation'] <- 0.06

  siteDat$a_sat[siteDat$type == 'Sloped without Accumulation'] <- 0.35
  siteDat$a_fc[siteDat$type == 'Sloped without Accumulation'] <- 0.2
  siteDat$a_pwp[siteDat$type == 'Sloped without Accumulation'] <- 0.06
  siteDat$o_sat[siteDat$type == 'Sloped without Accumulation'] <- 0.35
  siteDat$o_fc[siteDat$type == 'Sloped without Accumulation'] <- 0.2
  siteDat$o_pwp[siteDat$type == 'Sloped without Accumulation'] <- 0.06

  return(siteDat)

}


## Extract site data
siteDat <- get_siteDat(locs)

write.csv(siteDat, "scripts/inputs/KLC/intermediate_data/Landsat_init_soildat_output.csv", quote = F, row.names = F)
