#
# #########################
# Purpose: Functions for Treeline Migration Study
# Author: Adrianna C. Foster
# Date: June, 2022
# R version 4.2.0 (2022-04-22) 'Vigorous Calisthenics'
# #########################
# #########################

make_site_grid <- function(outline, size_km, seed){
  # Create a site grid given and outline and grid size
  # Inputs:
  # outline      - shapefile of study area (SpatialPolygonsDataFrame)
  # size_km      - grid size (km)
  # seed         - random seed for making site IDs (integer)
  # Output:
  # sites_latlon - data frame with site IDs and latitude/longitude
  
  #project outline to Albers so that units are in m
  Albers.crs <- CRS(paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 ",
                           "+lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 ",
                           "+datum=NAD83 +units=m +no_defs "))
  outline.albers <- spTransform(outline, CRS = Albers.crs)
  
  #get bounding box
  bbox <- bbox(outline.albers)
  
  ## create grid of sites and transform back into original CRS
  grid <- get_grid(bbox, size_km, Albers.crs, outline.albers)
  grid <- spTransform(grid, crs(outline))
  
  ## make data frame and rename columns
  sites_latlon <- as.data.frame(grid) %>%
    dplyr::select(y.1, x.1) %>% 
    rename('latitude' = 'y.1', 'longitude' = 'x.1')
  
  ## create random siteIDs
  set.seed(seed)
  IDs <- seq(10001, 999999, 1)
  siteIDs <- sample(IDs, size = nrow(sites_latlon), replace = FALSE)
  
  ## create site column
  sites_latlon$site <- siteIDs
  sites_latlon <- dplyr::select(sites_latlon, site, latitude, longitude)
  
  return(sites_latlon)
  
}

get_grid <- function(bbox, size_km, Albers.crs, outline){
  # create a grid of evenly spaced points from a bounding box
  # Inputs:
  # bbox     - bounding box
  # size_km  - grid size (km)
  # Output:
  # grid.pts - grid of points (SpatialPointsDataFrame)
  
  #create points for x and y from bounding box
  x <- seq(from = bbox[1, 1], to = bbox[1, 2], by = size_km*1000)
  y <- seq(from = bbox[2, 1], to = bbox[2, 2], by = size_km*1000)
  
  #create grid of pairs of coordinates
  xy <- expand.grid(x = x, y = y)
  
  #identify projection and create SPDF
  grid.pts <- SpatialPointsDataFrame(coords = xy, data = xy, 
                                     proj4string = Albers.crs)
  
  ## crop to study area
  grid.pts <- crop(grid.pts, outline)
  
  return(grid.pts)
}

make_raster_grid <- function(outline, size_km){
  # create a raster grid from an outline and grid size
  # Inputs:
  # outline  - outline (SptialPolygonsDataFrame)
  # size_km  - grid size (km)
  # Output:
  # rp       - raster grid (raster)
  
  ## project outline to Albers so that units are in m
  Albers.crs <- CRS(paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 ",
                           "+lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 ",
                           "+datum=NAD83 +units=m +no_defs "))
  outline.albers <- spTransform(outline, CRS = Albers.crs)
  
  ## get bounding box
  bbox <- bbox(outline.albers)
  
  ## create points for x and y from bounding box
  x <- seq(from = bbox[1, 1], to = bbox[1, 2], by = size_km*1000)
  y <- seq(from = bbox[2, 1], to = bbox[2, 2], by = size_km*1000)
  
  ##create grid of pairs of coordinates
  xy <- expand.grid(x = x, y = y)
  
  ##identify projection and create SPDF
  grid.pts <- SpatialPointsDataFrame(coords = xy, data = xy, 
                                     proj4string = Albers.crs)
  
  ## turn grid into polygon
  gridded(grid.pts) <- TRUE
  grid.pts$id <- row.names(grid.pts)
  
  ## create raster from grid
  ras <- raster(ext = extent(grid.pts), crs = Albers.crs, 
                res = c(size_km*1000, size_km*1000))
  
  ## rasterize the grid polygons
  rp <- rasterize(grid.pts, ras)
  
  return(rp[['id']])
  
  
}

get_i_val <- function(row, col, nrows){
  # calculate 1-D index
  # Inputs:
  # row      - row number
  # col      - column number
  # nrows    - total number of rows
  # Output:
  # ondDindex - 1-D index value
  
  oneDindex <- ((row-1) * nrows) + col
  
  return(oneDindex)
}

get_row_col_ind <- function(sites_latlon, ras_grid){
  # create columns for row, column, and 1-D index
  # Inputs:
  # sites_latlon - data frame with site IDs, latitude, and longitude
  # rp           - raster grid with same grid size as sites_latlon
  # Output:
  # siteDat     - data frame with site IDs, latitude, longitude, row, col, and
  #               1-D index
  
  ## create SPDF from input grid of lat/lons
  slocs <- SpatialPointsDataFrame(coords = sites_latlon[, c('longitude','latitude')], 
                                  data = sites_latlon, 
                                  proj4string = CRS("+init=epsg:4326"))
  
  ## transform slocs into the raster grid CRS
  slocs <- spTransform(slocs, crs(ras_grid))
  
  ## make raster of rows and columns
  cells <- seq(1, ncell(ras_grid))
  
  row_ras <- ras_grid
  row_ras[] <- rowFromCell(ras_grid, cells)
  
  col_ras <- ras_grid
  col_ras[] <- colFromCell(ras_grid, cells)
  
  ## extract row/cell values
  rows <- extract(row_ras, slocs)
  cols <- extract(col_ras, slocs)
  slocs@data <- cbind(slocs@data, rows, cols)
  
  ## turn into data frame
  siteDat <- as.data.frame(slocs)
  
  ## calculate 1-D index value
  nrows <- max(siteDat$rows)
  siteDat$ind <- get_i_val(siteDat$rows, siteDat$cols, nrows)
  
  ## check for duplicates, otherwise return siteDat
  ndupes <- sum(duplicated(siteDat$ind))
  if (ndupes > 0){
    print('WARNING: DUPLICATED INDICES')
    return(NULL)
  } else {
    return(dplyr::select(siteDat, site:cols, ind))
  }
  
}

capwords <- function(s, strict = FALSE) {
  # capitalize just first letter in each word of a character
  # Inputs:
  # s      - input character
  # strict - uncapitalize every other letter? 
  # Output:
  # out    - character with first letter capitalized
  
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  out <- sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
  return(out)
}

get_grid_ids <- function(dat, grid){
  
  #convert to spatial points df and transform polygon to correct crs
  coordinates(dat) <- ~longitude+latitude
  proj4string(dat) <- CRS("+init=epsg:4326")
  grid_t <- spTransform(grid, crs(dat))

  #get polygon IDs for each site
  over.data <- over(dat, grid_t[,'id'])
  dat@data$id <- over.data$id
  
  #covert back to dataframe and get rid of latitude and longitude
  sp_df <- as.data.frame(dat) %>%
    dplyr::select(site, latitude, longitude, id)
  
  return(sp_df)
  
}

make_spatial <- function(dat, sp_df, grid){
  
  dat <- merge(dat, dplyr::select(sp_df, site, id), by = 'site')
  
  #subset polygon to only areas where we have overlapping sites
  ids <- unique(dat$id)
  grid_sub <- subset(grid, grid@data$id %in% ids)
  
  #merge data - DO NOT use grid_sub@data, this will corrupt the spatial data!!!
  grid_sub <- merge(grid_sub, dat, by = 'id', all.x = FALSE)
  
  return(grid_sub)
  
}

make_species_spatial <- function(dat, sp_df, grid, yr, val) {
  
  dat_cast <- filter(dat, year == yr) %>%
    reshape2::dcast(tag + siteID + year ~ species,
                    value.var = val) %>%
    dplyr::rename('site' = 'siteID')
  
  grid_sub <- make_spatial(dat_cast, sp_df, grid)
  
  return(grid_sub)

}

param_site <- function(ex_sites){
  # parameterize a UVAFME2018_site.csv file from an input set of extracted
  # site and soil values
  # Inputs:
  # ex_sites - data frame with extracted site values
  # Output:
  # site_dat - data frame that can be used as an input UVAFME2018_site.csv file
  
  ## first grab relevant data from extracted values
  site_dat <- ex_sites %>%
    filter(!is.na(region) & !is.na(ex_obd) & !is.na(slope) & !is.na(ex_twi)) %>%
    mutate(a_fc = ex_afc/100,
           a_pwp = ex_apwp/100,
           o_fc = ex_ofc/100,
           o_pwp = ex_opwp/100,
           a_sat = a_pwp + a_fc,
           o_sat = o_pwp + o_fc,
           a_bd = ex_abd*10,
           o_bd = ex_obd*10,
           itxt = ifelse(ex_sand > 60, 0, 
                         ifelse(ex_sand > 40, 1, 2))) %>%
    mutate(hum_init = 55.3, A_depth = 1, wprob = 0, stand_age = 100,
           gcm_year = 100, flow = 0, 
           region = gsub(' ', '_', capwords(region, strict = T))) %>%
    dplyr::rename('row' = 'rows', 'col' = 'cols')
  
  ## apply some rules
  site_dat <- site_dat %>%
    mutate(norm_twi = ex_twi/132,
           type = ifelse(slope <= 2.0 & norm_twi > 0.75, 'flat_accum',
                         ifelse(slope <= 2.0 & norm_twi <= 0.75, 'flat_no_accum',
                                ifelse(slope > 2.0 & norm_twi > 0.75, 
                                       'sloped_accum', 'sloped_no_accum')))) %>%
    mutate(a_sat = ifelse(type == 'sloped_no_accum', 0.35, a_sat),
           a_fc = ifelse(type == 'sloped_no_accum', 0.2, a_fc),
           a_pwp = ifelse(type == 'sloped_no_accum', 0.06, a_pwp),
           o_sat = ifelse(type == 'sloped_no_accum', 0.35, o_sat),
           o_fc = ifelse(type == 'sloped_no_accum', 0.2, o_fc),
           o_pwp = ifelse(type == 'sloped_no_accum', 0.038, o_fc),
           itxt = ifelse(type == 'sloped_no_accum', 1, itxt),
           hum_init = ifelse(type == 'sloped_no_accum', 55.3/2, hum_init),
           
           
           a_sat = ifelse(type %in%  c('flat_no_accum', 'sloped_accum'), 0.44, a_sat),
           a_fc = ifelse(type %in% c('flat_no_accum', 'sloped_accum'), 0.29, a_fc),
           a_pwp = ifelse(type %in% c('flat_no_accum', 'sloped_accum'), 0.06, a_pwp),
           o_sat = ifelse(type %in% c('flat_no_accum', 'sloped_accum'), 0.38, o_sat),
           o_fc = ifelse(type %in% c('flat_no_accum', 'sloped_accum'), 0.3, o_fc),
           o_pwp = ifelse(type %in% c('flat_no_accum', 'sloped_accum'), 0.038, o_pwp),
           itxt = ifelse(type %in% c('flat_no_accum', 'sloped_accum'), 1, itxt),
           hum_init = ifelse(type %in%  c('flat_no_accum', 'sloped_accum'), 55.3*0.75, 
                             hum_init),
           
           a_sat = ifelse(type == 'flat_accum', 0.53, a_sat),
           a_fc = ifelse(type == 'flat_accum', 0.38, a_fc),
           a_pwp = ifelse(type == 'flat_accum', 0.06, a_pwp),
           o_sat = ifelse(type == 'flat_accum', 0.53, o_sat),
           o_fc = ifelse(type == 'flat_accum', 0.38, o_fc),
           o_pwp = ifelse(type == 'flat_accum', 0.038, o_pwp),
           itxt = ifelse(type == 'flat_accum', 2, itxt)) %>%
    dplyr::select(site, latitude, longitude, name, region, elevation, slope, 
                  aspect, a_sat, a_fc, a_pwp, o_sat, o_fc, o_pwp,
                  a_bd, o_bd, itxt, hum_init, A_depth, wprob,
                  stand_age, gcm_year, flow, row, col, ind)
  
  return(site_dat)
}

param_clim <- function(climNA, cldmn_dir, cldsd_dir, windmn_dir){
  # parameterize UVAFME climate files from input values
  # Inputs:
  # climNA     - output from Climate NA software for sites
  # cldmn_dir  - directory with mean monthly cloudiness rasters
  # cldsd_dir  - directory with sd monthly cloudiness rasters
  # windmn_dir - directory with mean monthly wind speed rasters
  # Outputs:
  # out        - list of data frames with input climate files for UVAFME
  
  
  ## get rid of columns we don't want
  climNA <- dplyr::select(climNA, Year:Longitude, Tmax01:Tmin12, PPT01:PPT12, 
                         RH01:RH12)
  
  ## set -9999 to NA
  climNA[climNA == -9999] <- NA
  
  ## group by id1 and get mean temp and precip
  climNA_means <- group_by(climNA, id1) %>%
    summarize(tmin_jan = mean(Tmin01, na.rm = TRUE), 
              tmin_feb = mean(Tmin02, na.rm = TRUE),
              tmin_mar = mean(Tmin03, na.rm = TRUE), 
              tmin_apr = mean(Tmin04, na.rm = TRUE),
              tmin_may = mean(Tmin05, na.rm = TRUE),
              tmin_jun = mean(Tmin06, na.rm = TRUE),
              tmin_jul = mean(Tmin07, na.rm = TRUE),
              tmin_aug = mean(Tmin08, na.rm = TRUE),
              tmin_sep = mean(Tmin09, na.rm = TRUE),
              tmin_oct = mean(Tmin10, na.rm = TRUE),
              tmin_nov = mean(Tmin11, na.rm = TRUE),
              tmin_dec = mean(Tmin12, na.rm = TRUE),
              
              tmax_jan = mean(Tmax01, na.rm = TRUE),
              tmax_feb = mean(Tmax02, na.rm = TRUE),
              tmax_mar = mean(Tmax03, na.rm = TRUE),
              tmax_apr = mean(Tmax04, na.rm = TRUE),
              tmax_may = mean(Tmax05, na.rm = TRUE),
              tmax_jun = mean(Tmax06, na.rm = TRUE),
              tmax_jul = mean(Tmax07, na.rm = TRUE),
              tmax_aug = mean(Tmax08, na.rm = TRUE),
              tmax_sep = mean(Tmax09, na.rm = TRUE),
              tmax_oct = mean(Tmax10, na.rm = TRUE),
              tmax_nov = mean(Tmax11, na.rm = TRUE),
              tmax_dec = mean(Tmax12, na.rm = TRUE),
              
              prcp_jan = mean(PPT01, na.rm = TRUE),
              prcp_feb = mean(PPT02, na.rm = TRUE),
              prcp_mar = mean(PPT03, na.rm = TRUE),
              prcp_apr = mean(PPT04, na.rm = TRUE),
              prcp_may = mean(PPT05, na.rm = TRUE),
              prcp_jun = mean(PPT06, na.rm = TRUE),
              prcp_jul = mean(PPT07, na.rm = TRUE),
              prcp_aug = mean(PPT08, na.rm = TRUE),
              prcp_sep = mean(PPT09, na.rm = TRUE),
              prcp_oct = mean(PPT10, na.rm = TRUE),
              prcp_nov = mean(PPT11, na.rm = TRUE),
              prcp_dec = mean(PPT12, na.rm = TRUE))
  
  ## group by id1 and get mean RH
  climNA_RHmeans <- group_by(climNA, id1) %>%
    summarize(rh_jan = mean(RH01, na.rm = TRUE),
              rh_feb = mean(RH02, na.rm = TRUE),
              rh_mar = mean(RH03, na.rm = TRUE),
              rh_apr = mean(RH04, na.rm = TRUE),
              rh_may = mean(RH05, na.rm = TRUE),
              rh_jun = mean(RH06, na.rm = TRUE),
              rh_jul = mean(RH07, na.rm = TRUE),
              rh_aug = mean(RH08, na.rm = TRUE),
              rh_sep = mean(RH09, na.rm = TRUE),
              rh_oct = mean(RH10, na.rm = TRUE),
              rh_nov = mean(RH11, na.rm = TRUE),
              rh_dec = mean(RH12, na.rm = TRUE))

  
  ## group by id1 and get sd of temp and precip
  climNA_sds <-  group_by(climNA, id1) %>%
    summarize(tmin_std_jan = sd(Tmin01, na.rm = TRUE),
              tmin_std_feb = sd(Tmin02, na.rm = TRUE),
              tmin_std_mar = sd(Tmin03, na.rm = TRUE),
              tmin_std_apr = sd(Tmin04, na.rm = TRUE),
              tmin_std_may = sd(Tmin05, na.rm = TRUE), 
              tmin_std_jun = sd(Tmin06, na.rm = TRUE), 
              tmin_std_jul = sd(Tmin07, na.rm = TRUE), 
              tmin_std_aug = sd(Tmin08, na.rm = TRUE),
              tmin_std_sep = sd(Tmin09, na.rm = TRUE),
              tmin_std_oct = sd(Tmin10, na.rm = TRUE),
              tmin_std_nov = sd(Tmin11, na.rm = TRUE),
              tmin_std_dec = sd(Tmin12, na.rm = TRUE),
                          
              tmax_std_jan = sd(Tmax01, na.rm = TRUE), 
              tmax_std_feb = sd(Tmax02, na.rm = TRUE),
              tmax_std_mar = sd(Tmax03, na.rm = TRUE), 
              tmax_std_apr = sd(Tmax04, na.rm = TRUE),
              tmax_std_may = sd(Tmax05, na.rm = TRUE), 
              tmax_std_jun = sd(Tmax06, na.rm = TRUE), 
              tmax_std_jul = sd(Tmax07, na.rm = TRUE), 
              tmax_std_aug = sd(Tmax08, na.rm = TRUE),
              tmax_std_sep = sd(Tmax09, na.rm = TRUE),
              tmax_std_oct = sd(Tmax10, na.rm = TRUE),
              tmax_std_nov = sd(Tmax11, na.rm = TRUE),
              tmax_std_dec = sd(Tmax12, na.rm = TRUE),
                          
              prcp_std_jan = sd(PPT01, na.rm = TRUE), 
              prcp_std_feb = sd(PPT02, na.rm = TRUE),
              prcp_std_mar = sd(PPT03, na.rm = TRUE), 
              prcp_std_apr = sd(PPT04, na.rm = TRUE),
              prcp_std_may = sd(PPT05, na.rm = TRUE), 
              prcp_std_jun = sd(PPT06, na.rm = TRUE), 
              prcp_std_jul = sd(PPT07, na.rm = TRUE), 
              prcp_std_aug = sd(PPT08, na.rm = TRUE),
              prcp_std_sep = sd(PPT09, na.rm = TRUE),
              prcp_std_oct = sd(PPT10, na.rm = TRUE),
              prcp_std_nov = sd(PPT11, na.rm = TRUE),
              prcp_std_dec = sd(PPT12, na.rm = TRUE))
  
  ## group by id1 and get sd of RH
  climNA_RHsds <-  group_by(climNA, id1) %>%
    summarize(rh_std_jan = sd(RH01, na.rm = TRUE),
              rh_std_feb = sd(RH02, na.rm = TRUE),
              rh_std_mar = sd(RH03, na.rm = TRUE),
              rh_std_apr = sd(RH04, na.rm = TRUE),
              rh_std_may = sd(RH05, na.rm = TRUE),
              rh_std_jun = sd(RH06, na.rm = TRUE),
              rh_std_jul = sd(RH07, na.rm = TRUE),
              rh_std_aug = sd(RH08, na.rm = TRUE),
              rh_std_sep = sd(RH09, na.rm = TRUE),
              rh_std_oct = sd(RH10, na.rm = TRUE),
              rh_std_nov = sd(RH11, na.rm = TRUE),
              rh_std_dec = sd(RH12, na.rm = TRUE))

  ## get site lat and long
  site_locs <- dplyr::select(climNA, id1, Latitude, Longitude) %>% 
    distinct_all()
  
  ## merge climate data with site locations to get latitude/longitude
  climNA_means <- merge(climNA_means, site_locs, by = 'id1')
  climNA_sds <- merge(climNA_sds, site_locs, by = 'id1')
  
  ## rename columns and reorder
  climNA_means <- rename(climNA_means, 'site' = 'id1', 'latitude' = 'Latitude',
                         'longitude' = 'Longitude') %>% 
    dplyr::select(site, latitude, longitude, tmin_jan:prcp_dec)
  
  climNA_sds <- rename(climNA_sds, 'site' = 'id1', 'latitude' = 'Latitude',
                       'longitude' = 'Longitude') %>% 
    dplyr::select(site, latitude, longitude, tmin_std_jan:prcp_std_dec)
  
  ## Monthly cloudiness input---
  
  ## get input cloud tif files
  cldmns <- list.files(cldmn_dir, pattern = '.tif$')
  cldsds <- list.files(cldsd_dir, pattern = '.tif$')
  
  ## read in raster stack of cloud means and sds
  cldmn <- stack()
  for (r in 1:length(cldmns)){
    cld_tmp <- raster(paste(cldmn_dir, cldmns[r], sep = '/'))
    cldmn <- stack(cldmn, cld_tmp)
  }
  
  cldsd <- stack()
  for (r in 1:length(cldsds)){
    cld_tmp <- raster(paste(cldsd_dir, cldsds[r], sep = '/'))
    cldsd <- stack(cldsd, cld_tmp)
  }
  
  
  ## convert site info to spatial points dataframe
  siteLocs_sp <- SpatialPointsDataFrame(coords = site_locs[,c(3,2)], 
                                        data = site_locs, 
                                        proj4string = crs(cldmn))
  
  ## extract cloud information at site locations
  ex_mn <- extract(cldmn, siteLocs_sp)
  ex_sd <- extract(cldsd, siteLocs_sp)
  siteLocs_sp@data <- cbind(siteLocs_sp@data, data.frame(ex_mn), 
                            data.frame(ex_sd))
  
  ## convert back to dataframe and split into mean and std columns
  cld_mn <- data.frame(siteLocs_sp) %>% dplyr::select(id1:cld_12mean)
  cld_sd <- data.frame(siteLocs_sp) %>% dplyr::select(id1:Longitude, 
                                                      cld_01std:cld_12std)
  
  ## rename columns
  cld_mn <- rename(cld_mn, 'site' = 'id1', 'latitude' = 'Latitude', 
                   'longitude' = 'Longitude', 'cld_jan' = 'cld_01mean',
                   'cld_feb' = 'cld_02mean', 'cld_mar' = 'cld_03mean',
                   'cld_apr' = 'cld_04mean', 'cld_may' = 'cld_05mean',
                   'cld_jun' = 'cld_06mean', 'cld_jul' = 'cld_07mean',
                   'cld_aug' = 'cld_08mean', 'cld_sep' = 'cld_09mean',
                   'cld_oct' = 'cld_10mean', 'cld_nov' = 'cld_11mean',
                   'cld_dec' = 'cld_12mean')
  
  cld_sd <- rename(cld_sd, 'site' = 'id1', 'latitude' = 'Latitude', 
                   'longitude' = 'Longitude', 'cld_jan_std' = 'cld_01std',
                   'cld_feb_std' = 'cld_02std', 'cld_mar_std' = 'cld_03std',
                   'cld_apr_std' = 'cld_04std', 'cld_may_std' = 'cld_05std',
                   'cld_jun_std' = 'cld_06std', 'cld_jul_std' = 'cld_07std',
                   'cld_aug_std' = 'cld_08std', 'cld_sep_std' = 'cld_09std',
                   'cld_oct_std' = 'cld_10std', 'cld_nov_std' = 'cld_11std',
                   'cld_dec_std' = 'cld_12std')
  
  ## Monthly wind input---
  
  windmns <- list.files(windmn_dir, pattern = '.tif$')
  
  ## read in stack of wind means 
  windmn <- stack()
  for (r in 1:length(windmns)){
    wnd_tmp <- raster(paste(windmn_dir, windmns[r], sep = '/'))
    windmn <- stack(windmn, wnd_tmp)
  }
  
  siteLocs_sp_t <- spTransform(siteLocs_sp, crs(windmn))
  
  ex_wind <- extract(windmn, siteLocs_sp_t)
  siteLocs_sp_t@data <- cbind(siteLocs_sp_t@data, data.frame(ex_wind))
  wind_mn <- data.frame(siteLocs_sp_t) %>% 
    dplyr::select(id1,wc2.1_30s_wind_01:wc2.1_30s_wind_12)
  
  ## merge cloud, RH, and wind data
  
  extra_mn <- merge(cld_mn, climNA_RHmeans, by.x = 'site', by.y = 'id1') %>%
    merge(wind_mn, by.x = 'site', by.y = 'id1')
  
  ## rename columns
  extra_mn <- rename(extra_mn, 'wind_jan' = 'wc2.1_30s_wind_01',
                     'wind_feb' = 'wc2.1_30s_wind_02', 
                     'wind_mar' = 'wc2.1_30s_wind_03',
                     'wind_apr' = 'wc2.1_30s_wind_04', 
                     'wind_may' = 'wc2.1_30s_wind_05',
                     'wind_jun' = 'wc2.1_30s_wind_06', 
                     'wind_jul' = 'wc2.1_30s_wind_07',
                     'wind_aug' = 'wc2.1_30s_wind_08', 
                     'wind_sep' = 'wc2.1_30s_wind_09',
                     'wind_oct' = 'wc2.1_30s_wind_10', 
                     'wind_nov' = 'wc2.1_30s_wind_11',
                     'wind_dec' = 'wc2.1_30s_wind_12')
  
  extra_sd <- merge(cld_sd, climNA_RHsds, by.x = 'site', by.y = 'id1')
  
  extra_mn <- dplyr::select(extra_mn, site:longitude, cld_jan:wind_dec)
  extra_sd <- dplyr::select(extra_sd, site:longitude, cld_jan_std:rh_std_dec)
  
  clim_means_check <- reshape2::melt(climNA_means, id.vars = 'site',
                                     measure.vars = c(4:39),
                                     variable.name = 'var_month') %>%
    mutate(var = sapply(strsplit(as.character(var_month), '_'), `[`, 1),
           month = sapply(strsplit(as.character(var_month), '_'), `[`, 2))
  clim_sds_check <- reshape2::melt(climNA_sds, id.vars = 'site',
                                     measure.vars = c(4:39),
                                     variable.name = 'var_month') %>%
    mutate(var = sapply(strsplit(as.character(var_month), '_'), `[`, 1),
           month = sapply(strsplit(as.character(var_month), '_'), `[`, 2))
  
  clim_ex_means_check <- reshape2::melt(extra_mn, id.vars = 'site',
                                     measure.vars = c(4:39),
                                     variable.name = 'var_month') %>%
    mutate(var = sapply(strsplit(as.character(var_month), '_'), `[`, 1),
           month = sapply(strsplit(as.character(var_month), '_'), `[`, 2))
  clim_ex_sds_check <- reshape2::melt(extra_sd, id.vars = 'site',
                                   measure.vars = c(4:27),
                                   variable.name = 'var_month') %>%
    mutate(var = sapply(strsplit(as.character(var_month), '_'), `[`, 1),
           month = sapply(strsplit(as.character(var_month), '_'), `[`, 2))
  
  temp_mean_bad <- filter(clim_means_check, var %in% c('tmin', 'tmax')) %>%
    filter(is.na(value))
  precip_mean_bad <- filter(clim_means_check, var == 'prcp') %>%
    filter(is.na(value) | value < 0)
  temp_sd_bad <- filter(clim_sds_check, var %in% c('tmin', 'tmax')) %>%
    filter(is.na(value) | value < 0)
  precip_sd_bad <- filter(clim_sds_check, var == 'prcp') %>%
    filter(is.na(value) | value < 0)
  
  ex_mean_bad <- filter(clim_ex_means_check, var %in% c('cld', 'rh', 'wind')) %>%
    filter(is.na(value) | value < 0)
  
  ex_sd_bad <- filter(clim_ex_means_check, var %in% c('cld', 'rh')) %>%
    filter(is.na(value) | value < 0)
  
  bad_sites <- unique(c(temp_mean_bad$site, precip_mean_bad$site,
                        temp_sd_bad$site, precip_sd_bad$site,
                        ex_mean_bad$site, ex_sd_bad$site))
  

  climNA_means <- filter(climNA_means, !(site %in% bad_sites))
  climNA_sds <- filter(climNA_sds, !(site %in% bad_sites))
  extra_mn <- filter(extra_mn, !(site %in% bad_sites))
  extra_sd <- filter(extra_sd, !(site %in% bad_sites))
  
  
  out <- list(clim_mn = climNA_means, clim_sd = climNA_sds, 
              extra_mn = extra_mn, extra_sd = extra_sd,
              bad_sites = bad_sites)
  
  return(out)

}

param_lightning <- function(lightning_dir, lightning_file, locs){
  # parameterize UVAFME lightning file from input values
  # Inputs:
  # lightning_dir   - directory with mean monthly lightning strike density
  # lightning_file  - shape file with mean monthly lightning strike density
  # locs            - site locations data frame
  # Outputs:
  # lightning_dat   - data frame that can serve as input lightning file for 
  #                   UVAFME
  
  shp <- readOGR(dsn = lightning_dir, layer = lightning_file)
  
  #convert site locations to a spatial points dataframe
  siteLocs_sp <- SpatialPointsDataFrame(coords = 
                                          locs[,c('longitude','latitude')], 
                                        data = locs,
                                        proj4string = CRS("+init=epsg:4326"))
  
  #extract site location data and convert NAs (i.e. data not present) to 0
  ex <- over(spTransform(siteLocs_sp, crs(shp)), shp)
  
  #cbind with original data frame
  siteLocs_sp@data <- cbind(siteLocs_sp@data, ex)
  
  lightning_dat <- as.data.frame(siteLocs_sp)
  
  lightning_dat <- dplyr::select(lightning_dat, site:longitude, 
                                 strmn_jan:strsd_dec)
  
  
  ## get rid of unwanted columns, set NA values to 0
  lightning_dat$na_count <- apply(lightning_dat, 1, function(x) sum(is.na(x)))
  
  lightning_dat <- filter(lightning_dat, na_count < 24) %>%
    dplyr::select(site:strsd_dec)
  lightning_dat[is.na(lightning_dat)] <- 0
  
  return(lightning_dat)
  
}

param_sitelist <- function(locs){
  # parameterize UVAFME sitelist file from input values
  # Inputs:
  # locs     - site locations data frame
  # Outputs:
  # sitelist - data frame that can serve as input sitelist file for UVAFME
  
  sitelist <- dplyr::select(locs, site)
  sitelist$runID <- row.names(sitelist)
  sitelist$altitude <- NA
  sitelist$management <- NA
  sitelist$year_management <- NA
  sitelist$shearblading <- NA
  sitelist$prescribed_burn <- NA
  sitelist$thinning <- NA
  sitelist$thin_perc <- NA
  sitelist$sel_cutting <- NA
  sitelist$rotation_time <- NA
  sitelist$dbh_thresh <- NA
  sitelist$fire_year <- NA
  sitelist$fire_day <- NA
  sitelist$fire_wind <- NA
  sitelist$ffmc <- NA
  sitelist$dmc <- NA
  sitelist$end <- NA
  
  return(sitelist)
  
}

param_rangelist <- function(locs, specRangeDir, specs_to_include, 
                            cavm_dir, cavm_treeline){
  ## function to parameterize rangelist from input species range maps 
  ## Inputs:
  ##  locs            - site locations, should have site, latitude, longitude
  ##  cavm_dir        - directory with CAVM treeline
  ##  cavm_treeline   - shapefile for CAVM treeline
  ##  specRangeDir    - directory with shapefiles of rangemaps
  ## specs_to_include - list of species to add to rangelist
  
  ## get directory and species names of available range maps
  dirs <- list.files(specRangeDir, full.names = TRUE)
  specs <- paste0(toupper(substr(basename(dirs), 1, 4)),
                  substr(basename(dirs), 5, 8))
  
  ## convert to spdf
  wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  locs_spdf <- SpatialPointsDataFrame(coords = locs[,c('longitude','latitude')],
                                      data = locs, 
                                      proj4string = wgs)
  
  cavm <- readOGR(cavm_dir, cavm_treeline)

  ## loop through and pull species range information from the maps
  for (s in 1:length(specs_to_include)){
    
    ## get current colnames
    cnames <- colnames(locs_spdf@data)
    
    ## read in specific range map and extract
    id <- which(specs == specs_to_include[s])
    
    ## if it doesn't exist, skip
    if (length(id) > 0){
      map <- readOGR(dsn = dirs[id], stringsAsFactors = FALSE)
      if (is.na(crs(map))){
        proj4string(map) <- CRS("+init=epsg:4267")
      }
      
      ex_spec <- as.numeric(over(spTransform(locs_spdf, crs(map)), map)$CODE)
      
      ## add to data
      locs_spdf@data <- cbind(locs_spdf@data, ex_spec)
      
      ## update colnames
      colnames(locs_spdf@data) <- c(cnames, specs_to_include[s])
    } else {
      locs_spdf@data[,specs_to_include[s]] <- NA
    }
  }
  
  ex_cavm <- as.numeric(over(spTransform(locs_spdf, crs(cavm)), 
                             cavm)$WATER)
  locs_spdf@data <- cbind(locs_spdf@data, ex_cavm)

  range_map <- as.data.frame(locs_spdf) %>%
    dplyr::select(-longitude.1, -latitude.1)
  
  ## add species that don't have a range map but we know are present
  range_map$DASIfrut <- 1
  range_map$LEDUgroe <- 1
  
  range_map[is.na(range_map)] <- 0
  
  range_map$BETUneoa[range_map$ex_cavm != 2] <- 0
  range_map$LARIlari[range_map$ex_cavm != 2] <- 0
  range_map$PICEglau[range_map$ex_cavm != 2] <- 0
  range_map$PICEmari[range_map$ex_cavm != 2] <- 0
  range_map$POPUbals[range_map$ex_cavm != 2] <- 0
  range_map$POPUtrem[range_map$ex_cavm != 2] <- 0
  
  spec_sums <- colSums(range_map[,c(4:ncol(range_map))])
  spec_include <- names(spec_sums[spec_sums > 0])
  
  range_map <- dplyr::select(range_map, all_of(c('site', 'latitude',
                                                 'longitude', spec_include))) %>%
    dplyr::select(-ex_cavm)
  
  return(range_map)
  
  
}

make_inputdirs_from_list <- function(maindir = NA, indir_main, input_prefix, 
                                     sitelists){
  #function to create input directories for UVAFME batch simulations for
  #a specific folder code
  #maindir - main I/O directory
  #indir_main - main Input directory
  #input_prefix - input folder code prefix
  #sitelists - list of sitelists
  
  #get name of main I/O directory (if using) and main input directory and 
  #create if they don't exist
  if (!is.na(maindir)){
    if (!(file.exists(maindir))){
      dir.create(maindir)
    }
    if (!file.exists(file.path(maindir, indir_main))){
      dir.create(file.path(maindir, indir_main))
    }
  } else {
    if (!file.exists(file.path(indir_main))){
      dir.create(file.path(indir_main))
    }
  }
  
  for (i in 1:length(sitelists)){
    slist_code <- gsub('UVAFME2018_sitelist_|.csv', '', basename(sitelists[i]))
    if (!is.na(maindir)){
      dir.create(file.path(maindir, indir_main, paste0(input_prefix, '_', 
                                                       slist_code)))
    } else {
      dir.create(file.path(indir_main, paste0(input_prefix, '_', slist_code)))
    }
  }
}

make_outputdirs_from_list <- function(maindir = NA, outdir_main, output_prefix, 
                                     sitelists){
  #function to create output directories for UVAFME batch simulations for
  #a specific folder code
  #maindir - main I/O directory
  #outdir_main - main output directory
  #output_prefix - output folder code prefix
  #sitelists - list of sitelists
  
  #get name of main I/O directory (if using) and main input directory and 
  #create if they don't exist
  if (!is.na(maindir)){
    if (!(file.exists(maindir))){
      dir.create(maindir)
    }
    if (!file.exists(file.path(maindir, outdir_main))){
      dir.create(file.path(maindir, outdir_main))
    }
  } else {
    if (!file.exists(file.path(outdir_main))){
      dir.create(file.path(outdir_main))
    }
  }
  
  for (i in 1:length(sitelists)){
    slist_code <- gsub('UVAFME2018_sitelist_|.csv', '', basename(sitelists[i]))
    if (!is.na(maindir)){
      dir.create(file.path(maindir, outdir_main, paste0(output_prefix, '_', 
                                                       slist_code)))
    } else {
      dir.create(file.path(outdir_main, paste0(output_prefix, '_', slist_code)))
    }
  }
}

make_fileLists_from_list <- function(maindir = NA, flist_dir, flist_name, inputmain, 
                           indir_main, outdir_main, input_prefix, output_prefix,
                           sitelists, wdir){
  #function to create file lists for UVAFME batch simulations for a
  #flist_dir - directory for filelists
  #flist_name - name of file list
  #inputmain - main input directory with main files
  #input_prefix - input folder prefix
  #output_prefix - output folder prefix
  #sitelists - list of sitelists
  #wdir - working directory
  
  if (!is.na(maindir)){
    if (!(file.exists(maindir))){
      dir.create(maindir)
    }
    if (!(file.exists(file.path(maindir, flist_dir)))){
      dir.create(file.path(maindir, flist_dir))
    }
  } else {
    if (!(file.exists(file.path(flist_dir)))){
      dir.create(file.path(flist_dir))
    }
  }
  
  if (!is.na(maindir)){
    wdir <- file.path(wdir, maindir)
  }
  
  #create and write file lists
  l1 <- '&filenames'
  l2 <- paste0("input_directory='", file.path(wdir, indir_main, inputmain), "'")
  #l3 #output directory
  l4 <- paste0("climate_directory='", file.path(wdir, indir_main, inputmain), "'")
  l5 <- paste0("site_directory='", file.path(wdir, indir_main, inputmain), "'")
  #l6 #sitelist directory
  #l7 #GCM directory
  l8 <- paste0("rt_directory='", file.path(wdir, indir_main, inputmain), "'")
  l9 <- paste0("speclist_directory='", file.path(wdir, indir_main, inputmain), "'")
  l10 <- '/'
  
  #for each sitelist
  for (i in 1:length(sitelists)){
    
     slist_code <- gsub('UVAFME2018_sitelist_|.csv', '', 
                        basename(sitelists[i]))
     
    #open file to write file list to
    if (!is.na(maindir)){
      fname <- paste0(maindir, '/', flist_dir, '/', flist_name, '_', 
                      slist_code, '.txt')
    } else {
      fname <- paste0(flist_dir, '/', flist_name, '_', slist_code, '.txt')
    }

    fileConn <- file(fname)
    
    #create lines
    l3 <- paste0("output_directory='", 
                 paste0(wdir, '/', outdir_main, '/', 
                           paste(output_prefix, slist_code, sep = '_'), "'"))
    l6 <- paste0("sitelist_directory='", 
                 paste0(wdir, '/', indir_main, '/',
                           paste(input_prefix, slist_code, sep = '_'), "'"))
    l7 <- paste0("GCM_directory='", 
                 paste0(wdir, '/', indir_main, '/',
                           paste(input_prefix, slist_code, sep = '_'), "'"))
    
    flines <- c(l1, l2, l3, l4, l5, l6, l7, l8, l9, l10)
    writeLines(flines, fileConn)
    close(fileConn)
  }
  
  
}

distribute_files_from_list <- function(maindir = NA, gcm_main = NA, 
                                       indir_main, input_prefix, sitelists){
  #function to distribute sitelist and potentially GCM files from a list
  #maindir - main directory
  #gcm_main - main gcm file with all sites
  #indir_main - main output folder
  #input_prefix - input folder prefix
  #sitelists - list of site lists
  
  for (i in 1:length(sitelists)){
    slist_code <- gsub('UVAFME2018_sitelist_|.csv', '', basename(sitelists[i]))
    if (!is.na(maindir)){
      indir <- file.path(maindir, indir_main, paste0(input_prefix, '_', 
                                                        slist_code))
    } else {
      indir <- file.path(indir_main, paste0(input_prefix, '_', slist_code))
    }
    
    #subset sitelist and gcm file
    sitelist_sub <- read.csv(sitelists[i])
    
    if (!is.na(gcm_main)){
      gcmfile_sub <- subset(gcm_main, gcm_main$site %in% 
                              unique(sitelist_sub$site))
    }
    
    #write to folder
    write.csv(sitelist_sub, paste0(indir, '/', 'UVAFME2018_sitelist.csv'), 
              row.names = FALSE, na = '')
    
    if (!is.na(gcm_main)){
      write.csv(gcmfile_sub, paste0(indir, '/', 
                                    'UVAFME2018_climate_GCM.csv'), 
                row.names = FALSE, na = '')
    }
    
  }
  
}

make_script_from_list <- function(maindir = NA, flist_dir, flist_name, jcode,
                                   wdir, exec, jobscript_dir, sitelists,
                                  list_file){
  #script to create shell script for submitting UVAFME batch jobs to SLURM
  #Inputs:
  #maindir: main directory
  #flist_dir: directory where filelists are located (chr)
  #flist_name: filelist or nodelist name prefix (chr)
  #jcode: job code name
  #wdir: working directory (chr)
  #exec: executable name (chr)
  #jobscript_dir: name of jobscript folder (chr)
  #sitelists - list of sitelists
  
  #create job name and file name
  job_file <- paste0('job_', jcode, '_${code}.slurm')
  job_name <- paste0(jcode, '_${code}')
  flist <- paste0(flist_dir, '/', flist_name, '_${code}.txt')
  
  jobs <- length(sitelists)
  
  #create and write lines of shell script
  l1 <- '#!/bin/bash'
  l2 <- '#'
  l3 <- ' '
  l4 <- 'arr=()'
  l5 <- 'while IFS= read -r line; do'
  l6 <- 'arr+=("$line")'
  l7 <- paste0('done<', list_file)
  
  l8 <- paste0('for i in ${arr[@]}')
  l9 <- 'do'
  l10 <- paste0('code=$i')
  
  l11 <- paste0('echo "#!/bin/bash" >> ', job_file)
  l12 <- paste0('echo "#SBATCH --job-name=', job_name, '" >> ', job_file)
  l13 <- paste0('echo "#SBATCH -o outjob_', job_name, '" >> ', job_file)
  l14 <- paste0('echo "#SBATCH -e erjob_', job_name, '" >> ', job_file)
  l15 <- paste0('echo "#SBATCH --workdir=', wdir, '" >> ', job_file)
  l16 <- paste0('echo "', exec, ' ', flist, '" >> ', job_file)
  l17 <- paste0('sbatch ', job_file)
  l18 <- 'done'
  
  lines <- c(l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14,
             l15, l16, l17, l18)
  
  #create folder for job scripts
  if (!is.na(maindir)){
    if (!dir.exists(maindir)){
      dir.create(maindir)
    }
    if(!dir.exists(file.path(maindir, jobscript_dir))){
      dir.create(file.path(maindir, jobscript_dir))
    }
  } else {
    if(!dir.exists(file.path(jobscript_dir))){
      dir.create(jobscript_dir)
    }
  }
  
  #create file to write to
  output_code <- jcode
  if (!is.na(maindir)){
    fname <- file.path(maindir, 
                       paste0(jobscript_dir, '/runjob_', output_code, '.sh'))
  } else {
    fname <- paste0(jobscript_dir, '/runjob_', output_code, '.sh')
  }
  
  fileConn <- file(fname)
  
  writeLines(lines, fileConn)
  close(fileConn)
  
  if (!is.na(maindir)){
    fileConn <- file(paste0(maindir, '/', jobscript_dir, '/', list_file))
  } else {
    fileConn <- file(paste0(jobscript_dir, '/', list_file))
  }

  codes <- gsub('UVAFME2018_sitelist_|.csv', '', basename(sitelists))
  writeLines(codes, fileConn)
  close(fileConn)
}

compile_outdat <- function(out.dir, out_prefix = 'output_data_', fname) {
  #compile output data across multiple folders
  #Inputs:
  #out.dir    - main output directory
  #out_prefix - prefix for sub-directories
  #fname      - name of file to compile
  #Output:
  #out.dat    - compiled data
  
  ## main output files
  out.dirs <- list.files(out.dir, full.names = T)
  
  out.dat <- data.frame()
  for (i in 1:length(out.dirs)){
    file.tag <- gsub(out_prefix, '', basename(out.dirs[i]))
    dat <- fread(paste0(out.dirs[i], '/', fname))
    dat[dat==-999] <- NA
    dat$tag <- file.tag
    out.dat <- rbind(out.dat, dat)
  }
  
  return(out.dat)
  
}



graph_map <- function(dat, maxV, year) {
  
  ## create color palette 
  pal <- colorRampPalette(brewer.pal(9, 'Greens'))
  scf <- scale_fill_gradientn(colours = pal(100), limits = c(0, maxV))
  
  dat.fort <- fortify(dat, region = 'site')
  dat.fort <- merge(dat.fort, dat, by.x = 'id', by.y = 'site')
  dat.fort$PICEglau[dat.fort$PICEglau <= 0.01] <- NA
  
  g <- ggplot(dat.fort, aes(x = long, y = lat, group = id)) +
    geom_polygon(aes(fill = PICEglau)) + 
    scf +
    facet_wrap(.~tag, scales = 'free') +
    theme_bw() +
    theme(axis.text = element_blank()) +
    ggtitle(paste0('Year ', year))
  
  return(g)
  
}
