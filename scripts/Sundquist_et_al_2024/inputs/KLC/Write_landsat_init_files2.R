
library(terra)
library(dplyr)
library(reshape2)
library(data.table)
library(scales)
library(corrplot)
library(ggplot2)
# library('raster')
# library('rgdal')
# library('gdalUtilities')
# library('gdalUtils')
# library('sp')
# library('foreign')
# source("https://raw.githubusercontent.com/gonzalezivan90/SDG15_indicators/main/Ecuador_fragmentacion/R_03_tabuleRaster.R")
source("scripts/analyses/KLC/Utilities.R")

LandsatDir = "raw/LandsatTreeSpecies/"


# Read in data
age = rast(paste0(LandsatDir, "forestage_2019/CA_forest_age_2019.tif"))  #years
BA = rast(paste0(LandsatDir, "basalarea_2015/CA_forest_basal_area_2015_NN.tif")) #m2/ha
CC = rast(paste0(LandsatDir, "canopycover_2015/CA_forest_percentage_first_returns_above_2m_2015.tif"))
Domspec = rast(paste0(LandsatDir, "lead_species_2019/CA_forest_lead_tree_species.tif"))

# Crop and reproject data
oldcrs = crs(age)
wgs84   = "+proj=longlat +datum=WGS84"
tile = ext(-128,-127,59,60) #WGS84 coords
tempext = ext(-1.85e6, -1.65e6,1.5e6, 1.8e6) # Preliminary crop to make projection faster

age = crop(age, tempext)
BA = crop(BA, tempext)
CC = crop(CC, tempext)
Domspec = crop(Domspec, tempext)

# Summarize data on 1km grid
nad83_14n = "+proj=aea +lat_0=45 +lon_0=-126 +lat_1=50 +lat_2=58.5 +x_0=1000000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"
tile2 = project(tile, from = wgs84, to =nad83_14n)
r = rast(ext(Domspec), res = 1000, crs = oldcrs) 
r = rast(ext = tile2, res = 1000, crs = nad83_14n)

Domspec = project(Domspec, nad83_14n, method = "near")
sfir = Domspec
sfir[sfir != 3,] <- 0
sfirsum = resample(sfir, r, "sum")
sfircount = sfirsum/3

sfircount = project(sfircount, wgs84, method = "near")
sfircount = as.points(sfircount)
sfircount = as.data.frame(sfircount, geom = "XY")
colnames(sfircount) = c("sfir","x","y")
sfircount$x = round(sfircount$x, 4)
sfircount$y = round(sfircount$y, 4)

wsp = Domspec
wsp[wsp != 17,] <- 0
wspsum = resample(wsp, r, "sum")
wspcount = wspsum/17

wspcount = project(wspcount, wgs84, method = "near")
wspcount = as.points(wspcount)
wspcount = as.data.frame(wspcount, geom = "XY")
colnames(wspcount) = c("wsp","x","y")
wspcount$x = round(wspcount$x, 4)
wspcount$y = round(wspcount$y, 4)

bsp = Domspec
bsp[bsp != 18,] <- 0
bspsum = resample(bsp, r, "sum")
bspcount = bspsum/18

bspcount = project(bspcount, wgs84, method = "near")
bspcount = as.points(bspcount)
bspcount = as.data.frame(bspcount, geom = "XY")
colnames(bspcount) = c("bsp","x","y")
bspcount$x = round(bspcount$x, 4)
bspcount$y = round(bspcount$y, 4)

lpine = Domspec
lpine[lpine != 23,] <- 0
lpinesum = resample(lpine, r, "sum")
lpinecount = lpinesum/23

lpinecount = project(lpinecount, wgs84, method = "near")
lpinecount = as.points(lpinecount)
lpinecount = as.data.frame(lpinecount, geom = "XY")
colnames(lpinecount) = c("lpine","x","y")
lpinecount$x = round(lpinecount$x, 4)
lpinecount$y = round(lpinecount$y, 4)

aspen = Domspec
aspen[aspen != 29,] <- 0
aspensum = resample(aspen, r, "sum")
aspencount = aspensum/29

aspencount = project(aspencount, wgs84, method = "near")
aspencount = as.points(aspencount)
aspencount = as.data.frame(aspencount, geom = "XY")
colnames(aspencount) = c("aspen","x","y")
aspencount$x = round(aspencount$x, 4)
aspencount$y = round(aspencount$y, 4)


specprops = merge(sfircount, wspcount, by = c("x","y"))
specprops = merge(specprops, bspcount, by = c("x","y"))
specprops = merge(specprops, lpinecount, by = c("x","y"))
specprops = merge(specprops, aspencount, by = c("x","y"))

specprops$sum = specprops$sfir + specprops$wsp + specprops$bsp + specprops$lpine + specprops$aspen 
ggplot(specprops, aes(x, y, color = sum)) + geom_point()
# check: sum should be less than 1111 (30m pixels to 1km2)... looks exactly right, yay!

# Trim the edges
specprops = specprops[specprops$x < -127 & specprops$x > -128 & specprops$y > 59 & specprops$y < 60,]

specprops[,3:7] = specprops[,3:7]/specprops[,8]

# Visualize
temp = melt(specprops[,1:7], id.vars = c("x","y"))
temp$variable = factor(temp$variable, levels = c("sfir","wsp","bsp","lpine","aspen"),
                          labels = c("Abies lasiocarpa","Picea glauca","Picea mariana","Pinus contorta","Populous tremuloides"))
g<- ggplot(temp, aes(x, y, color = value)) + facet_wrap(~variable) + geom_point() + scale_color_continuous(type = "viridis")
ggsave(paste0(figDir, "species_props.png"), g,  width = 8, height = 6)

# Add basal area, age, canopy cover to the data frame. Using near to project because interpolation might produce unrealistic combinations of species, BA, CC, etc. 
BA[BA == 0,] <- NA
BA = project(BA, nad83_14n, method = "near")
BAmean = resample(BA, r)

BAmean = project(BAmean, wgs84, method = "near")
BAmean = as.points(BAmean)
BAmean = as.data.frame(BAmean, geom = "XY")
colnames(BAmean) = c("BA","x","y")
BAmean$x = round(BAmean$x, 4)
BAmean$y = round(BAmean$y, 4)

# age[age == 0,] <- NA
age = project(age, nad83_14n, method = "near")
agemean = resample(age, r)

agemean = project(agemean, wgs84, method = "near")
agemean = as.points(agemean)
agemean = as.data.frame(agemean, geom = "XY")
colnames(agemean) = c("age","x","y")
agemean$x = round(agemean$x, 4)
agemean$y = round(agemean$y, 4)

CC[CC == 0,] <- NA
CC = project(CC, nad83_14n, method = "near")
CCmean = resample(CC, r)

CCmean = project(CCmean, wgs84, method = "near")
CCmean = as.points(CCmean)
CCmean = as.data.frame(CCmean, geom = "XY")
colnames(CCmean) = c("CC","x","y")
CCmean$x = round(CCmean$x, 4)
CCmean$y = round(CCmean$y, 4)
  
dat = merge(specprops, BAmean, by = c("x","y"))
dat = merge(dat, agemean, by = c("x","y"))
dat = merge(dat, CCmean, by = c("x","y"))

# Visualize
temp = dat
temp$BA = temp$BA/max(temp$BA)
temp$age = temp$age/max(temp$age)
temp$CC = temp$CC/max(temp$CC)
temp = melt(temp[,c(1:7, 9:11)], id.vars = c("x","y"))
ggplot(temp, aes(x, y, color = value)) + facet_wrap(~variable) + geom_point()


# Classify dat into stand phases
dat$phase = ""
dat[dat$age < 10 | dat$CC < 50,]$phase <- "SI"
dat[dat$CC >= 50 & dat$age <= 34 & dat$phase != "SI",]$phase <- "SE"
dat[dat$age > 34 & dat$aspen >= 0.5,]$phase <- "SE"
dat[dat$age > 34 & dat$aspen < 0.5,]$phase <- "CT" # aspen is equal to deciduous fraction
dat[dat$age >= 100,]$phase <- "GD"


# Begin to init_sp input file
init_sp = dat[,c(1:7, 9)]
init_sp$id = 1:nrow(init_sp)
colnames(init_sp) = c("lon","lat","ABIElasi","PICEglau","PICEmari","PINUcont","POPUtrem","BA","siteID")

# Begin init_dbh input file
set.seed(111)
init_dbh = init_sp[,c(1,2,9)]
SSDparams = fread("scripts/inputs/KLC/intermediate_data/SSDparams.csv")
cats = c("less5","5to10","10to20","20to30","30to40","40to50","50to60","60to70","70to80","80up")
sp = colnames(init_sp)[3:7]
cols = paste0(sp[1], "_", cats)
cols = c(cols, paste0(sp[2], "_", cats))
cols = c(cols, paste0(sp[3], "_", cats))
cols = c(cols, paste0(sp[4], "_", cats))
cols = c(cols, paste0(sp[5], "_", cats))

# warning: if dat or init_sp get re-sorted this will be fucked up 

for(i in 3:7) { #columns of init_sp corresponding with species proportions ... except ABIES, which isn't in params yet
  spec = colnames(init_sp)[i]
  specbins = data.frame(matrix(NA, ncol=13, nrow=nrow(dat), dimnames=list(NULL, c("siteID","lon","lat",  paste0(spec, "_", cats)))))

  for(j in 1:nrow(dat)) {
    params = SSDparams[SSDparams$species == spec & SSDparams$phase == dat[j,]$phase,]
    samp = rweibull(10000, params$shape, params$scale)
    bins = samp %>% cut(breaks = c(0,5,10,20,30,40,50,60,70,80,2000))
    specbins[j,] = c(init_sp[j,c("siteID","lon","lat")], table(bins))
  }
  init_dbh = merge(init_dbh, specbins, by= c("siteID","lon","lat"))
}

fwrite(init_sp, "scripts/inputs/KLC/intermediate_data/landsat_init_sp_test.csv")
fwrite(init_dbh, "scripts/inputs/KLC/intermediate_data/landsat_init_dbh_test.csv")
