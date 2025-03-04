
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(segmented)
library(stringr)
library(EnvStats)

outDir = "output_data/KLC/SSDtest/"
LandsatDir = "raw/LandsatTreeSpecies/"

specdat = fread(paste0(outDir, "Species_Data.csv"))
specdat[specdat == -999] <- NA
specdat = specdat[,-c(28:33)]
specdat = na.omit(specdat)
SSD = melt(specdat[,c(1,3,5:16)], id.vars = c("species","siteID","year"))

dead = fread(paste0(outDir, "Dead_Species_Data.csv"))
dead[dead == -999] <- NA
dead = na.omit(dead)

plots = fread(paste0(outDir, "Total_Plot_Values.csv"))
# Fixing bug in output... this has been fixed in fortran code to prevent future errors
plots[,21:24] = plots[,21:24] - plots[,13:16]

mgmt <- fread(paste0(outDir, "Active_Management_Data.csv"))

mgmt2 = data.frame(matrix(ncol = 4, nrow = 1))
colnames(mgmt2) = c("siteID","classification","start","end")
ix = 1
for(i in unique(mgmt$siteID)) {
  y  = min(mgmt[mgmt$siteID == i,]$year)
  mgmt2 = rbind(mgmt2, c(i, mgmt[mgmt$siteID == i & mgmt$year == y,]$classification, y, NA))
  ix = ix + 1
  for (j in (y + 1):max(mgmt[mgmt$siteID == i,]$year)) {
    if(mgmt[mgmt$siteID == i & mgmt$year == j,]$classification != mgmt[mgmt$siteID == i & mgmt$year == (j-1),]$classification){
      mgmt2[ix,]$end = j - 1
      mgmt2 = rbind(mgmt2, c(i, mgmt[mgmt$siteID == i & mgmt$year == j,]$classification, j, NA))
      ix = ix + 1
    } 
    if (j == max(mgmt$year)) {
      mgmt2[ix,]$end = max(mgmt$year)
    }
  }
}

mgmt2 = mgmt2[-1,]
mgmt2$start = as.double(mgmt2$start); mgmt2$end = as.double(mgmt2$end)


site = 1

#Classifications ~ biomass
p1 <- ggplot(subset(specdat, siteID %in% site)) +
  geom_area(aes(fill = species, x = year, y = total_biomC)) + 
  geom_vline(data = subset(mgmt2, siteID %in% site), aes(xintercept = end)) +
  geom_text(data = subset(mgmt2, siteID %in% site), aes(label = classification, x = (start + end)/2, y = 60))



p2 <-ggplot(subset(SSD, year == 20 & siteID %in% site), aes(x = variable, y = value, fill = species)) + geom_col() + 
  facet_wrap(~species)

ggarrange(p1, p2, nrow = 2)

p3 <- ggplot(subset(dead, siteID %in% site), aes(x = year, y = shade_death, color = species)) + geom_line()
p3 <- ggplot(subset(specdat, siteID %in% site), aes(x = year, y = shade_resp, color = species)) + geom_line()

plots_shade = melt(plots[,c(1,3,21:24)], id.vars = c("year","siteID"))
p3 <- ggplot(subset(plots_shade, siteID %in% site), aes(x = year, y = value, color = variable)) + geom_line()

p4 <- ggplot(subset(plots, siteID %in% site), aes(x = year, y = LAI_1)) + geom_line()


# Deciduous fraction by basal area
specdat$decid = 0
specdat[specdat$species %in% c("ALNUsuni","ALNUtenu","POPUbals","POPUtrem","BETUneoa","SALIscou")]$decid = 1
DF = specdat %>% group_by(year, siteID, decid) %>% summarize(sum(basal_area))

DF = merge(DF, plots[,c("siteID","year","basal_area")], by = c("siteID","year"))
DF = DF[DF$decid == 1,]
DF$DF = DF$`sum(basal_area)`/DF$basal_area

p5 <- ggplot(subset(DF, siteID %in% site), aes(x = year, y = DF)) + geom_line()

p6 <- ggplot(subset(dead, siteID %in% site), aes(x = year, y = total_biomC, color = species)) + geom_line()
p6 <- ggplot(subset(dead, siteID %in% site), aes(x = year, y = mean_diam, color = species)) + geom_line()

p6 <- ggplot(subset(specdat, siteID %in% site), aes(x = year, y = qmean_diam, color = species)) + geom_line()

# Looking for canopy transition
QDBH = specdat %>% group_by(year, siteID, decid) %>% summarize(sum(qmean_diam))
colnames(QDBH) = c("year","siteID","decid","qdbh")
p7 <- ggplot(subset(QDBH, siteID %in% site), aes(year, QDBH, color = as.factor(decid))) + geom_line()

# Run through all anticipated types of sites
site = 1
site = 101
site = 201
site = 301
#site = 401; 
site = 402 #401 is wacky

print(site)
p1 <- ggplot(subset(specdat, siteID %in% site)) +
  geom_area(aes(fill = species, x = year, y = total_biomC)) + 
  geom_vline(data = subset(mgmt2, siteID %in% site), aes(xintercept = end)) +
  geom_text(data = subset(mgmt2, siteID %in% site), aes(label = classification, x = (start + end)/2, y = 60))
p4 <- ggplot(subset(plots, siteID %in% site), aes(x = year, y = LAI_1)) + geom_line()

p7 <- ggplot(subset(QDBH, siteID %in% site), aes(year, QDBH, color = as.factor(decid))) + geom_line()
ggarrange(p1, p4, p7)

# Adding stand phase variable to specdat
# 1 = stand initiation, 2 = stem exclusion, 3 = canopy transition, 4 = gap dynamics 

site = 1
SI = 1 # start of stand initialization

lm.temp = lm(LAI_1 ~ year, data = subset(plots, plots$siteID %in% site))
bp = mgmt2[mgmt2$siteID == 1,]$end
bp = bp[-length(bp)]
bp = c(10,50)
lm.seg = segmented(lm.temp, ~year, psi = bp)
plot.segmented(lm.seg)
points(plots[plots$siteID == site,]$year, plots[plots$siteID == site,]$LAI_1)

SE = as.integer(lm.seg$psi[1,"Est."]) # start of stem exclusion

QDBH = dcast(as.data.frame(QDBH), siteID + year~ decid)
colnames(QDBH) = c("siteID","year","con","dec")
temp = subset(QDBH, year > SE & siteID == site)
CT = min(temp[temp$con > temp$dec,]$year) # start of canopy transition

spec = unique(specdat[specdat$siteID %in% site,]$species)
GD = NA
for(y in 1:130) {
  for(s in spec) {
    temp = specdat[specdat$siteID == site & specdat$species == s & specdat$year %in% y:(y+20),]
    lm.temp = lm(total_biomC ~ year, data = temp)
    if(abs(summary(lm.temp)$coefficients[2,1]) > 0.05) {
      break
    }
    # print(paste(y, s, lm.temp$coefficients))
    
    if(s == last(spec)) {
      GD = y # start of gap dynamics
      break
    }
  }
  if(!is.na(GD)) {
    break
  }
}

p8 <- p1 + geom_vline(xintercept = c(SI, SE, CT, GD), color = "red")

# Now turn this into function that will assign stand phase to every year site-by-site

assign_standphase = function(plotdat, specdat, qdbhdat) {
  SI = 151; SE = 151; CT = 151; GD = NA
  out = as.data.frame(cbind(1:150, rep(NA,150)))
  colnames(out) = c("year","phase")
  
  SI = 1 # start of stand initialization
  out[out$year >= SI,]$phase = "SI"
  
  lm.temp = lm(LAI_1 ~ year, data = plotdat)
  lm.seg = segmented(lm.temp, ~year, psi = c(10,50))

  SE = as.integer(lm.seg$psi[1,"Est."]) # start of stem exclusion
  out[out$year >= SE,]$phase = "SE"
  
  qdbhdat = subset(qdbhdat, year > SE)
  
  CT = min(qdbhdat[qdbhdat$con > qdbhdat$dec,]$year, 151) # start of canopy transition
  out[out$year >= CT,]$phase = "CT"
  
  spec = unique(specdat$species)
  GD = NA
  for(y in 1:130) {
    for(s in spec) {
      temp = specdat[specdat$species == s & specdat$year %in% y:(y+20),]
      lm.temp = lm(total_biomC ~ year, data = temp)
      if(abs(summary(lm.temp)$coefficients[2,1]) > 0.05) {
        break
      }
      # print(paste(y, s, lm.temp$coefficients))
      
      if(s == last(spec)) {
        GD = y # start of gap dynamics
        break
      }
    }
    if(!is.na(GD)) {
      out[out$year >= GD,]$phase = "GD"
      break
    }
  }
  
  
  return(out)
}

#Remove empty stands... input must have been wacky!
site = c(250 ,251, 309)
p1 <- ggplot(subset(specdat, siteID %in% site)) +
  geom_area(aes(fill = species, x = year, y = total_biomC)) + facet_wrap(~siteID)

specdat = specdat[!specdat$siteID %in% c(250, 251, 309),]

phases = data.frame(matrix(nrow = 0, ncol = 3))
colnames(phases) = c(c("siteID","year","phase"))
for(i in unique(specdat$siteID)) {
  plotdat = plots[plots$siteID == i,]
  specdat2 = specdat[specdat$siteID == i,]
  qdbhdat = QDBH[QDBH$siteID == i,]
  temp = assign_standphase(plotdat, specdat2, qdbhdat)
  temp$siteID = i
  phases = rbind(phases, temp)
  if(i %% 50 == 0) {
    print(".")
  }
}

temp = phases %>% group_by(phase, year) %>% summarize(count = n())

SSD = merge(SSD, phases, by = c("siteID","year"))
SSD = as.data.frame(SSD)
SSD$lower = as.integer((str_extract(SSD$variable, "^[0-9]{1,2}")))
SSD$upper = as.integer((str_extract(SSD$variable, "[0-9]{1,2}$")))
SSD$phase = factor(SSD$phase, levels = c("SI","SE","CT","GD"))

head(SSD)

#Doing runifs for each bin (e.g. 0-5, 5-10) and fit weibull off that
# 
# temp = subset(SSD, species == "POPUtrem" & siteID == 1 & year == 5)
expand <- function(df) {
  l = c()
  df = df[df$value > 0,]
  for(i in 1:nrow(df)) {
    l = c(l, runif(round(df[i,]$value), df[i,]$lower, df[i,]$upper))
  }
  return(l)
}

#Visualize lots of lines together
# Getting rid of shrubs
specdat$shape = -999; specdat$scale = -999
SSD = subset(SSD, !species %in% c("ALNUsinu","ALNUtenu","SALIscou"))
SSDparams = data.frame(matrix(nrow = 1, ncol = 6))
colnames(SSDparams) = c(c("phase","species","shape","scale","vshape","vscale"))
for(p in unique(SSD$phase)) {
  for(sp in unique(SSD$species)) {
    allcSSD = c()
    allshape = c()
    allscale = c()
    psp <- subset(SSD, phase == p & species == sp & value > 0.5)
    png(paste0("scripts/inputs/KLC/intermediate_data/inits_curves/",p,"_",sp,".png"), width = 960, height = 960)
    # size = nrow(psp %>% group_by(siteID, year))
    psp$plot = runif(nrow(psp), 1, 10)
    plot(1, type="n", xlab="", ylab="", xlim=c(0, 90), ylim=c(0, 0.5))
    for(s in unique(psp$siteID)) {
      psps = subset(psp, siteID == s) 
      for(y in unique(psps$year)) {
        df = subset(psps, year == y)
        cSSD = expand(df)
        # allcSSD = c(allcSSD, cSSD)
        if(length(cSSD) > 1 & df$plot[1] >= 9) {
          tempdist = eweibull(cSSD)
          curve(dweibull(x, shape=tempdist$parameters[1], scale = tempdist$parameters[2]), add = T, from=0, to=90, col = alpha("black", 0.2))
          allshape = c(allshape, tempdist$parameters[1])
          allscale = c(allscale, tempdist$parameters[2])
          specdat[specdat$species == sp & siteID == s & year == y,]$shape = tempdist$parameters[1]
          specdat[specdat$species == sp & siteID == s & year == y,]$scale = tempdist$parameters[2]
        }
      }
    }
    if(nrow(psp > 0)) {
      # pspdist = eweibull(allcSSD)
      # curve(dweibull(x, shape=pspdist$parameters[1], scale = pspdist$parameters[2]), add = T, from=0, to=90, col = "red", lwd = 3,  
      #       main = paste0(p, " ", sp))
      mshape = mean(allshape)
      mscale = mean(allscale)
      vshape <- var(allshape)
      vscale <- var(allscale)
      
      curve(dweibull(x, shape=mshape, scale = mscale), add = T, from=0, to=90, col = "red", lwd = 3, 
            main = paste0(p, " ", sp))
      SSDparams = rbind(SSDparams, c(p, sp, mshape, mscale, vshape, vscale))
    }
    dev.off()
  }
}

fwrite(SSDparams, "scripts/inputs/KLC/intermediate_data/SSDparams.csv")
fwrite(specdat, "scripts/inputs/KLC/intermediate_data/specdat_weibulls.csv")
fwrite(phases, "scripts/inputs/KLC/intermediate_data/phases.csv")


SSDparams = fread("scripts/inputs/KLC/intermediate_data/SSDparams.csv")

specdat = merge(specdat, phases, by = c("siteID","year"))
specdat[specdat == -999] <- NA
specdat = na.omit(specdat)
cor(specdat$basal_area, specdat$shape)
cor(specdat$basal_area, specdat$scale)


plots = fread("output_data/SSDtest/Total_Plot_Values.csv")
plots = merge(plots, phases, by = c("siteID","year"))
plots = merge(plots, DF, by = c("siteID","year"))

temp = plots %>% group_by(phase) %>% summarize(quantile(DF, na.rm = T))
plots$phase = factor(SSD$phase, levels = c("SI","SE","CT","GD"))

g<- ggplot(plots, aes(x = phase, y = DF, color = phase)) + geom_boxplot()
g<- ggplot(plots, aes(x = phase, y = year, color = phase)) + geom_boxplot()

#canopy transition curves by species
ggplot(data.frame(x = c(0,40)), aes(x)) + 
  stat_function(fun= dweibull, args = list("shape" = 1.29, "scale" = 4.85), color = "#FF7F00", linewidth = 2) +#lpine
  stat_function(fun= dweibull, args = list("shape" = 1.7, "scale" = 5.12), color = "#33A02C", linewidth = 2) + #birch
  stat_function(fun= dweibull, args = list("shape" = 1.12, "scale" = 5.74), color = "#FFFF99", linewidth = 2) + #aspen
  stat_function(fun= dweibull, args = list("shape" = 1.4, "scale" = 13.33), color = "#FB9A99", linewidth = 2) + #wsp
  stat_function(fun= dweibull, args = list("shape" = 1.17, "scale" = 6.76), color = "#1F78B4", linewidth = 2) + #bsp
  stat_function(fun= dweibull, args = list("shape" = 1.16, "scale" = 5.85), color = "#B2DF8A", linewidth = 2) + #balsam
  stat_function(fun= dweibull, args = list("shape" = 1.23, "scale" = 6.02), color = "#6A3D9A", linewidth = 2) + #fir
  theme_classic() + ylab("") + xlab("")

ggsave("scripts/inputs/KLC/intermediate_data/CT_species_weibulls.png", dpi = 600)






# import utilities for consistent species colors legend 
