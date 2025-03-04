### A script to explore Landsat-derived data products over Canada. These products will be converted to stand initialization data for UVAFME

library(terra)
library(dplyr)
library(reshape2)
library(data.table)
library(scales)
library(corrplot)

# LandsatDir = ""
# These data are very large and not provided in raw/; see raw/README.md for more information

# Read in data
age = rast(paste0(LandsatDir, "forestage_2019/CA_forest_age_2019.tif"))  #years
BA = rast(paste0(LandsatDir, "basalarea_2015/CA_forest_basal_area_2015_NN.tif")) #m2/ha
CC = rast(paste0(LandsatDir, "canopycover_2015/CA_forest_percentage_first_returns_above_2m_2015.tif"))
H = rast(paste0(LandsatDir, "forestheight_2015/CA_forest_elev_mean_2015_NN.tif"))
H_SD = rast(paste0(LandsatDir, "forestheight_sd_2015/CA_forest_elev_stddev_2015_NN.tif"))
Domspec = rast(paste0(LandsatDir, "lead_species_2019/CA_forest_lead_tree_species.tif"))
# Crop and reproject data



newcrs  = newcrs = "+proj=longlat +datum=WGS84"
tile = ext(-128,-127,59,60) #WGS84 coords
tempext = ext(-1.85e6, -1.65e6,1.5e6, 1.8e6) # Preliminary crop to make projection faster

age = crop(age, tempext)
age = project(age, newcrs, method = "near")
age = crop(age, tile)

BA = crop(BA, tempext)
BA = project(BA, newcrs)
BA = crop(BA, tile)

CC = crop(CC, tempext)
CC = project(CC, newcrs)
CC = crop(CC, tile)

H_SD = crop(H_SD, tempext)
H_SD = project(H_SD, newcrs)
H_SD = crop(H_SD, tile)

H = crop(H, tempext)
H = project(H, newcrs)
H = crop(H, tile)

Domspec = crop(Domspec, tempext)
Domspec = project(Domspec, newcrs, method = "near")
Domspec = crop(Domspec, tile)


# Convert to data frame
age = as.points(age)
age = as.data.frame(age, geom = "XY")
colnames(age) = c("age","x","y")
age$x = round(age$x, 4)
age$y = round(age$y, 4)

BA = as.points(BA)
BA = as.data.frame(BA, geom = "XY")
colnames(BA) = c("BA","x","y")
BA$x = round(BA$x, 4)
BA$y = round(BA$y, 4)

CC = as.points(CC)
CC = as.data.frame(CC, geom = "XY")
colnames(CC) = c("CC","x","y")
CC$x = round(CC$x, 4)
CC$y = round(CC$y, 4)

H_SD = as.points(H_SD)
H_SD = as.data.frame(H_SD, geom = "XY")
colnames(H_SD) = c("H_SD","x","y")
H_SD$x = round(H_SD$x, 4)
H_SD$y = round(H_SD$y, 4)

H = as.points(H)
H = as.data.frame(H, geom = "XY")
colnames(H) = c("H","x","y")
H$x = round(H$x, 4)
H$y = round(H$y, 4)

Domspec = as.points(Domspec)
Domspec = as.data.frame(Domspec, geom = "XY")
colnames(Domspec) = c("Domspec","x","y")
Domspec$x = round(Domspec$x, 4)
Domspec$y = round(Domspec$y, 4)

# Read in all species

specprobs = merge(BA, age, by = c("x","y"))
speclist = c("ABIE.AMA","ABIE.BAL","ABIE.LAS","ACER.MAC","ACER.RUB","ACER.SAH","ALNU.INC","ALNU.RUB",
             "BETU.ALL","BETU.PAP","CHAM.NOO","FRAX.NIG","LARI.LAR","LARI.OCC","PICE.ABI","PICE.ENG","PICE.GLA",
             "PICE.MAR","PICE.RUB","PICE.SIT","PINU.ALB","PINU.BAN","PINU.CON","PINU.PON","PINU.RES","PINU.STR",
             "POPU.BAL","POPU.GRA","POPU.TRE","PSEU.MEN","QUER.RUB","THUJ.OCC","THUJ.PLI","TSUG.CAN",
             "TSUG.HET","TSUG.MER","ULMU.AME")
for(s in speclist) {
  spec = rast(paste0(LandsatDir, "CA_forest_tree_species_probabilities/CA_forest_", s, "_2019.tif")) # percent cover
  spec = crop(spec, tempext)
  spec = project(spec, newcrs)
  spec = crop(spec, tile)
  spec = as.points(spec)
  spec = as.data.frame(spec, geom = "XY")
  colnames(spec) = c("prob","x","y")
  spec = spec[spec$prob > 0,]
  colnames(spec) = c(s,"x","y")
  if(nrow(spec) > 0) {
    spec$x = round(spec$x, 4)
    spec$y = round(spec$y, 4)
    specprobs = left_join(specprobs, spec, by = c("x","y"))
  }
}

specprobs = merge(specprobs, CC, by = c("x","y"))
specprobs = merge(specprobs, H_SD, by = c("x","y"))
specprobs = merge(specprobs, H, by = c("x","y"))
specprobs = merge(specprobs, Domspec, by = c("x","y"))

# specprobs = merge(specprobs, age, by = c("x","y"))
# specprobs = merge(specprobs, BA, by = c("x","y"))
# 
# # Save this!!! Took ages
fwrite(specprobs, paste0(LandsatDir, "tile/dat_59N_128W_specprob_age_BA.csv"))


dat = fread(paste0(LandsatDir, "tile/dat_59N_128W_specprob_age_BA.csv"))

dat <- dat %>% replace(is.na(.), 0)
corrplot(cor(dat))

# Kmeans clustering of data
library(cluster)
library(factoextra)
library(tidyverse)

dat2 <- scale(dat)

set.seed(123)

ix = sample(1:nrow(dat2), 10000)
dat2 = dat2[ix,3:11]


# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(dat2, k, iter.max = 20, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

k2 <- kmeans(dat2, centers = 2, nstart = 25)
k3 <- kmeans(dat2, centers = 3, nstart = 25)
k4 <- kmeans(dat2, centers = 4, nstart = 25)
k5 <- kmeans(dat2, centers = 5, nstart = 25)

# plots to compare
p1 <- fviz_cluster(k2, geom = "point", data = dat2) + ggtitle("k = 2")
p2 <- fviz_cluster(k3, geom = "point",  data = dat2) + ggtitle("k = 3")
p3 <- fviz_cluster(k4, geom = "point",  data = dat2) + ggtitle("k = 4")
p4 <- fviz_cluster(k5, geom = "point",  data = dat2) + ggtitle("k = 5")

library(gridExtra)
grid.arrange(p1, p2, p3, p4, nrow = 2)

dat2 = dat
dat2 = dat2[ix,]
dat2 = cbind(dat2, k5$cluster)

library(ggplot2)
ggplot(dat2, aes(x, y, color = as.factor(V2))) + geom_point()
dat2 %>% group_by(V2) %>% summarize_all(summary())


# Looking at CC as metric of stand phase


set.seed(123)

ix = sample(1:nrow(dat), 50000)
dat2 = dat[ix,3:12]

dat$age = as.integer(dat$age)
dat2 = dat %>% group_by(age) %>% sample_n(100, replace = F)
ss = smooth.spline(dat2$age, dat2$CC)
plot(dat2$age, dat2$CC) #,  col = alpha("black", 0.01))
lines(ss$x, ss$y, col = "red")
abline(v = 34, "blue")
abline(v = 54, "blue")

dat2 = dat %>% group_by(age) %>% sample_n(100, replace = F)
ss = smooth.spline(dat2$age, dat2$H_SD)
plot(dat2$age, dat2$H_SD) #,  col = alpha("black", 0.01))
lines(ss$x, ss$y, col = "red")


dat2 = dat %>% group_by(age) %>% sample_n(100, replace = F)
dat2 = na.omit(dat2)
ss = smooth.spline(dat2$age, dat2$NSD)
plot(dat2$age, dat2$NSD) #,  col = alpha("black", 0.01))
lines(ss$x, ss$y, col = "red")
