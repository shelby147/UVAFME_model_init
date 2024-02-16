#
# #########################
# Purpose: All figures and analyses for manuscript. Copying final versions for pub here so it's quick to update figs and stats if model is re-run
# Author: Shelby W. Sundquist
# Date: June, 2023
# R version 4.1.2
# #########################
# #########################
# Input format: txt, csv
# Output format: png, txt
# #########################
# 


library(data.table)
library(stringr)
library(sjmisc)
library(ggplot2)
library(survival)
library(ggsankey)
library(dplyr)
library(hash)
library(zoo)
library(glmnet)
library(survminer)
library(corrplot)
library(tidyr)
source("scripts/analyses/KLC/Utilities.R")
figDir = "scripts/analyses/KLC/figs/"

### Run CompileFiles.R 

f = "Active_Management_Data.csv"
mgmtdat = data.frame(matrix(ncol=length(cols[[f]]) + 1,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen"))))

for(c in clims) {
  path = paste0(outDir, c, "/", f)
    temp = fread(path)
  temp$climscen = c
  mgmtdat = rbind(mgmtdat, temp)
}

mgmt = mgmtdat

mgmtdat$classification = sapply(mgmtdat$classification, FUN = function(x) classes11hash[[toString(x)]])
mgmttemp = mgmtdat # snap of current dataset
years = seq(2015, 2100, 5)
mgmtdat$classification = factor(mgmtdat$classification)
mgmtdat =  as.data.frame(mgmtdat %>% group_by(siteID, climscen) %>% subset(year %in% years) %>% mutate(nextclass = lead(classification, n = 1)))
mgmtdat = as.data.frame(mgmtdat %>% group_by(siteID, climscen) %>% mutate(nextyear = lead(year, n = 1)))

mgmtdat = mgmt
mgmtdat$classification2 = ""
mgmtdat[mgmtdat$classification %in% c("BSP","WSPP","WSPS"),]$classification2 = "Spruce"
mgmtdat[mgmtdat$classification %in% c("ASP","ABI","BIR"),]$classification2 = "Deciduous"
mgmtdat[mgmtdat$classification %in% c("MIX","WSB","APIN"),]$classification2 = "Mixed forest"
mgmtdat[mgmtdat$classification %in% c("APIN"),]$classification2 = "Aspen/ pine"
mgmtdat[mgmtdat$classification %in% c("REP","OTH"),]$classification2 = "Reproductive/ other"
mgmtdat[mgmtdat$classification %in% c("LPIN"),]$classification2 = "Pine"

mgmtdat$classification = mgmtdat$classification2

mgmttemp = mgmtdat # snap of current dataset
years = c(2016, seq(2020, 2100, 5))
mgmtdat$classification = factor(mgmtdat$classification)
mgmtdat =  as.data.frame(mgmtdat %>% group_by(siteID, climscen) %>% subset(year %in% years) %>% mutate(nextclass = lead(classification, n = 1)))
mgmtdat = as.data.frame(mgmtdat %>% group_by(siteID, climscen) %>% mutate(nextyear = lead(year, n = 1)))


g <- ggplot(subset(mgmtdat), aes(x = year, next_x = nextyear,
                                 node = classification, next_node = nextclass,
                                 fill = factor(classification))) +
  facet_wrap(~climscen, nrow = 3) +
  geom_sankey(flow.alpha = 0.9, width = 0.5)  + theme_linedraw()+
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(colour = 'black'), axis.text = element_text(size = 15)) + 
  scale_fill_manual(breaks = classes6long, values = classpal6, name = "Forest classification") + 
  geom_vline(xintercept = c(2020,2040,2060,2080,2100), linetype = 3) +
  ylab("Number of sites") + xlab("Year") 
ggsave(paste0(figDir, "sankey_clims_broadclasses.png"), g, width = 10, height = 8)


g <- ggplot(subset(mgmtdat), aes(x = year, next_x = nextyear,
                                 node = classification, next_node = nextclass,
                                 fill = factor(classification))) +
  facet_wrap(~climscen, nrow = 3) +
  geom_sankey(flow.alpha = 0.9, width = 0.5)  + theme_linedraw()+
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(colour = 'black'), axis.text = element_text(size = 15)) + 
  scale_fill_manual(breaks = classes6long, values = classpal6, name = "Forest classification") + 
  geom_vline(xintercept = c(2020,2040,2060,2080,2100), linetype = 3) +
  ylab("Number of sites") + xlab("Year") +transparent
ggsave(paste0(figDir, "sankey_clims_broadclasses_transparent.png"), g, width = 10, height = 8)

### 2. 20-year landscape-wide biomass density and forest classifications between climate scenarios

f = "Active_Management_Data.csv"
mgmtdat = data.frame(matrix(ncol=length(cols[[f]]) + 1,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen"))))

for(c in clims) {
  path = paste0(outDir, c, "/", f)
  temp = fread(path)
  temp$climscen = c
  mgmtdat = rbind(mgmtdat, temp)
}

f = "Total_Plot_Values.csv"
plotsdat = data.frame(matrix(ncol=length(cols[[f]]) + 1,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen"))))
for(c in clims) {
  path = paste0(outDir, c, "/", f)
  temp = fread(path)
  temp$climscen = c
  plotsdat = rbind(plotsdat, temp)
}

dat = merge(plotsdat, mgmtdat, by = c("siteID","year","climscen"))
dat = as.data.frame(dat)
dat$classification2 = ""
dat[dat$classification %in% c("BSP","WSPP","WSPS"),]$classification2 = "Spruce"
dat[dat$classification %in% c("ASP","ABI","BIR"),]$classification2 = "Deciduous"
dat[dat$classification %in% c("MIX","WSB","APIN"),]$classification2 = "Mixed forest"
dat[dat$classification %in% c("APIN"),]$classification2 = "Aspen/ pine"
dat[dat$classification %in% c("REP","OTH"),]$classification2 = "Reproductive/ other"
dat[dat$classification %in% c("LPIN"),]$classification2 = "Pine"

dat$climscen2 = climnames(as.character(dat$climscen))
dat$climscen2 = factor(dat$climscen2, c("Historic climate","Moderate climate warming","Extreme climate warming"))

years = c(2020,2040,2060,2080,2100)

g <- ggplot(subset(dat, year %in% years)) + 
  geom_density(aes(x = total_biomC, fill = classification2, color = classification2,  alpha = 0.33), position = "stack") + 
  facet_grid(climscen2~year, switch = "y")  + ylim(0,0.25) + 
  scale_fill_manual(breaks = classes6long, values = classpal6, name = "Forest type") +
  scale_color_manual (breaks = classes6long, values = classpal6, name = "Forest type") + 
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        strip.text = element_text(colour = 'black')) +
  xlab("Biomass (tC/ha)") + ylab("Density") + scale_alpha(guide = 'none') + transparent

ggsave(paste0(figDir, "clim_biomass_class_20yrs.png"), g, width = 12, height = 8)


### 2.3 20-year relative changes in forest composition, biomass by forest type, and landscape-wide biomass
#Relative difference - baseline is historical climate scenario in corresponding year 
years = c(2016,2020,2040,2060,2080,2100)
foresttype_biomass_summary = matrix(nrow = 0, ncol = 6)
colnames(foresttype_biomass_summary) = c("classification", "mn_biomass","se_biomass","n", "year", "clim")

for(y in years) {
  for(c in clims) {
    ft_biomass = dat %>% group_by(classification2) %>% subset(year == y & climscen == c) %>% summarize(mean(total_biomC, na.rm = T))
    ft_biomass_se = dat %>% group_by(classification2) %>% subset(year == y & climscen == c) %>% summarize(sd(total_biomC, na.rm = T)/sqrt(n()))
    ft_sites = dat %>% group_by(classification2) %>% subset(year == y & climscen == c) %>% summarize(n())
    ft_biomass = cbind(ft_biomass, ft_biomass_se[,2])
    ft_biomass = cbind(ft_biomass, ft_sites[,2])
    colnames(ft_biomass) = c("classification", "mn_biomass","se_biomass","n")
    ft_biomass$clim = c
    ft_biomass$year = y
    foresttype_biomass_summary = rbind(foresttype_biomass_summary, ft_biomass)
  }
}

# fill in gaps
# foresttype_biomass_summary = rbind(foresttype_biomass_summary, c("Reproductive/ other", NA, NA, 0, "gcm45",2020))
foresttype_biomass_summary = rbind(foresttype_biomass_summary, list("Deciduous", NA, NA, 0, "gcm85",2100))
foresttype_biomass_summary = foresttype_biomass_summary[order(foresttype_biomass_summary$year, foresttype_biomass_summary$clim, foresttype_biomass_summary$classification),]

#biomass difference
foresttype_biomass_summary$diff = NA
diff = foresttype_biomass_summary[foresttype_biomass_summary$clim == "gcm45",]$mn_biomass - foresttype_biomass_summary[foresttype_biomass_summary$clim == "hist",]$mn_biomass
foresttype_biomass_summary[foresttype_biomass_summary$clim == "gcm45",]$diff = diff
diff = foresttype_biomass_summary[foresttype_biomass_summary$clim == "gcm85",]$mn_biomass - foresttype_biomass_summary[foresttype_biomass_summary$clim == "hist",]$mn_biomass
foresttype_biomass_summary[foresttype_biomass_summary$clim == "gcm85",]$diff = diff

foresttype_biomass_summary$rdiff = NA
rdiff = foresttype_biomass_summary[foresttype_biomass_summary$clim == "gcm45",]$mn_biomass / foresttype_biomass_summary[foresttype_biomass_summary$clim == "hist",]$mn_biomass
foresttype_biomass_summary[foresttype_biomass_summary$clim == "gcm45",]$rdiff = rdiff
rdiff = foresttype_biomass_summary[foresttype_biomass_summary$clim == "gcm85",]$mn_biomass / foresttype_biomass_summary[foresttype_biomass_summary$clim == "hist",]$mn_biomass
foresttype_biomass_summary[foresttype_biomass_summary$clim == "gcm85",]$rdiff = rdiff

foresttype_biomass_summary$rse_biomass = NA
rse_biomass = foresttype_biomass_summary[foresttype_biomass_summary$clim == "gcm45",]$se_biomass / foresttype_biomass_summary[foresttype_biomass_summary$clim == "hist",]$mn_biomass
foresttype_biomass_summary[foresttype_biomass_summary$clim == "gcm45",]$rse_biomass = rse_biomass
rse_biomass = foresttype_biomass_summary[foresttype_biomass_summary$clim == "gcm85",]$se_biomass / foresttype_biomass_summary[foresttype_biomass_summary$clim == "hist",]$mn_biomass
foresttype_biomass_summary[foresttype_biomass_summary$clim == "gcm85",]$rse_biomass = rse_biomass

#count difference
foresttype_biomass_summary$ndiff = NA
ndiff = foresttype_biomass_summary[foresttype_biomass_summary$clim == "gcm45",]$n - foresttype_biomass_summary[foresttype_biomass_summary$clim == "hist",]$n
foresttype_biomass_summary[foresttype_biomass_summary$clim == "gcm45",]$ndiff = ndiff
ndiff = foresttype_biomass_summary[foresttype_biomass_summary$clim == "gcm85",]$n - foresttype_biomass_summary[foresttype_biomass_summary$clim == "hist",]$n
foresttype_biomass_summary[foresttype_biomass_summary$clim == "gcm85",]$ndiff = ndiff

temp = foresttype_biomass_summary %>% select(classification, year, clim, n, mn_biomass, se_biomass)
temp$clim = factor(temp$clim, levels = c("hist","gcm45","gcm85"), labels = c("Historical climate","SSP2-4.5", "SSP5-8.5"))
colnames(temp) = c("Forest type","Year","Climate scenario","N","Biomass mn","Biomass se")

fwrite(temp, paste0(figDir, "supp_t5.csv"))

g<-ggplot(subset(foresttype_biomass_summary, clim != "hist" & n > 10 & classification != "Reproductive/ other" & year >= 2020), aes(x = year, y = rdiff, group = year,
                                                                                                                                    ymax = rdiff+ 1.96*rse_biomass, ymin = rdiff - 1.96*rse_biomass)) +
  geom_pointrange(aes(color = classification, fatten = n/1500, size = 1), position =position_dodge2(width = 7)) + geom_hline(yintercept = 1) +
  geom_point(shape = 21, aes(size = n/100), position = position_dodge2(width = 7)) +
  facet_wrap(~clim, nrow = 1) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(breaks = classes6long, values = classpal6, name = "Forest classification") +
  theme_linedraw() + ylab("Biomass change compared to historic climate") + xlab("Year") +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(colour = 'black'), 
        axis.text = element_text(size = 15)) + guides(size = "none",) 
ggsave(paste0(figDir, "/exogenous_classes_biomass_frequency_clims.png"), g, width = 12, height = 4)

g<-ggplot(subset(foresttype_biomass_summary, clim != "hist" & n > 10 & classification != "Reproductive/ other" & year >= 2020), aes(x = year, y = rdiff, group = year,
                                                                                                                                    ymax = rdiff+ 1.96*rse_biomass, ymin = rdiff - 1.96*rse_biomass)) +
  geom_pointrange(aes(color = classification, fatten = n/1500, size = 1), position =position_dodge2(width = 7)) + geom_hline(yintercept = 1) +
  geom_point(shape = 21, aes(size = n/100), position = position_dodge2(width = 7)) +
  facet_wrap(~clim, nrow = 1) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(breaks = classes6long, values = classpal6, name = "Forest classification") +
  theme_linedraw() + ylab("Biomass change compared to historic climate") + xlab("Year") +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(colour = 'black'), 
        axis.text = element_text(size = 15)) + guides(size = "none") + transparent
ggsave(paste0(figDir, "/exogenous_classes_biomass_frequency_clims_transparent.png"), g, width = 12, height = 4)

#Landscape-wide biomass
biomass_summary = matrix(nrow = 0, ncol = 4)
colnames(biomass_summary) = c("biomass","se_biomass","year", "clim")

years = c(2020,2060,2100)
for(y in years) {
  temp = (dat %>% group_by(climscen) %>% subset(year == y) %>% summarize(sum(total_biomC*100))) #convert ha to site area (25ha)
  temp2 = (dat %>% group_by(climscen) %>% subset(year == y) %>% summarize(sqrt(mean((total_biomC_sd*100)^2))/sqrt(n())))
  temp = cbind(temp, temp2[,2])
  colnames(temp) = c("clim","biomass","se_biomass")
  temp$year = y 
  biomass_summary = rbind(biomass_summary, temp)
}
biomass_summary$clim = factor(biomass_summary$clim, levels = c("hist","gcm45","gcm85"), labels = c("Historical climate","SSP2-4.5", "SSP5-8.5"))

ggplot(biomass_summary, aes(x = as.factor(year), y = biomass, fill = clim)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = biomass - 1.96*se_biomass, ymax = biomass + 1.96*se_biomass), position = position_dodge()) + 
  scale_fill_manual(breaks = c("Historical climate","SSP2-4.5", "SSP5-8.5"), values = c("#D1E5F0", "#F4A582", "#B2182B"), name = "Climate scenario") + 
  theme_linedraw() + theme(axis.text = element_text(size = 15))+ ylab("Canada tile-wide Biomass (tC)") + xlab("Year")

ggsave(paste0(figDir, "totalbiomass_clims.png"))

ggplot(biomass_summary, aes(x = as.factor(year), y = biomass, fill = clim)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = biomass - 1.96*se_biomass, ymax = biomass + 1.96*se_biomass), position = position_dodge()) + 
  scale_fill_manual(breaks = c("Historical climate","SSP2-4.5", "SSP5-8.5"), values = c("#D1E5F0", "#F4A582", "#B2182B"), name = "Climate scenario") + 
  theme_linedraw() + theme(axis.text = element_text(size = 15))+ ylab("Canada tile-wide Biomass (tC)") + xlab("Year") + transparent

ggsave(paste0(figDir, "totalbiomass_clims_transparent.png"))

rm(mgmtdat)
rm(plotsdat)
rm(firedat)
rm(dat)

# Spinup vs init figs
f = "Active_Management_Data.csv"
mgmtdat = data.frame(matrix(ncol=length(cols[[f]]) + 2,nrow=0, dimnames=list(NULL, c(cols[[f]], "iscen", "cscen"))))

for(c in clims) {
  path = paste0(outDir, "spinup/", c, "/", f)
  temp = fread(path)
  temp$cscen = c
  temp$iscen = "Spin-up"
  mgmtdat = rbind(mgmtdat, temp)
}
sites = unique(mgmtdat$siteID)

for(c in clims) {
  path = paste0(outDir, c, "/", f)
  temp = fread(path)
  temp = temp[temp$siteID %in% sites,]
  temp$cscen = c
  temp$iscen = "Initialization"
  mgmtdat = rbind(mgmtdat, temp)
}

f = "Total_Plot_Values.csv"
plotsdat = data.frame(matrix(ncol=length(cols[[f]]) + 2,nrow=0, dimnames=list(NULL, c(cols[[f]], "iscen", "cscen"))))

for(c in clims) {
  path = paste0(outDir, "spinup/", c, "/", f)
  temp = fread(path)
  temp$cscen = c
  temp$iscen = "Spin-up"
  plotsdat = rbind(plotsdat, temp)
}

for(c in clims) {
  path = paste0(outDir, c, "/", f)
  temp = fread(path)
  temp = temp[temp$siteID %in% sites,]
  temp$cscen = c
  temp$iscen = "Initialization"
  plotsdat = rbind(plotsdat, temp)
}

dat = merge(plotsdat, mgmtdat, by = c("siteID","year","iscen","cscen"))

mgmttemp = mgmtdat[mgmtdat$cscen == "hist",]
mgmttemp$classification = sapply(mgmttemp$classification, FUN = function(x) classes11hash[[toString(x)]])
years = seq(1913, 2100, 5)
mgmttemp$classification = factor(mgmttemp$classification)
mgmttemp$iscen = factor(mgmttemp$iscen, levels = c("Spin-up","Initialization"))
mgmttemp =  as.data.frame(mgmttemp %>% group_by(siteID, iscen) %>% subset(year %in% years) %>% mutate(nextclass = lead(classification, n = 1)))
mgmttemp = as.data.frame(mgmttemp %>% group_by(siteID,  iscen) %>% mutate(nextyear = lead(year, n = 1)))

g <- ggplot(subset(mgmttemp), aes(x = year, next_x = nextyear,
                                 node = classification, next_node = nextclass,
                                 fill = factor(classification))) +
  geom_sankey(flow.alpha = 0.9, width = 0.5) +
  facet_wrap(~iscen, nrow = 2) +
  theme_linedraw() + scale_fill_manual(breaks = classes11long, values = classpal11, name = "Forest classification") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(colour = 'black', size = 15),
        # legend.position = c(0.2,0.25),
        axis.text = element_text(size = 15)) +
  ylab("Number of sites") + xlab("Year") + ggtitle("UVAFME initialization vs spinup under historic climate")
ggsave(paste0(figDir, "init_spinup_hist_sankey.png"), g, width = 10, height = 8)

years = c(2020,2040, 2060,2080,2100)
endo_bm_diff = as.data.frame(dat %>% group_by(siteID, cscen, year) %>% subset(dat$year %in% years) %>% reframe(diff = diff(total_biomC)))
#diff is init - spinup; what is the exogenous bias from spinup?
endo_bm_diff$cscen = factor(endo_bm_diff$cscen, levels = c("hist","gcm45","gcm85"), labels = c("Historical climate","SSP2-4.5", "SSP5-8.5"))

ggplot(subset(endo_bm_diff, year %in% c(2020, 2060, 2100)), aes(x = as.factor(year), y = diff, fill = cscen)) +
  geom_boxplot(position = "dodge") + geom_hline(yintercept = 0) + 
  scale_fill_manual(breaks = c("Historical climate","SSP2-4.5", "SSP5-8.5"), values = c("#D1E5F0", "#F4A582", "#B2182B"), name = "Climate scenario") + 
  theme_linedraw() + theme(axis.text = element_text(size = 15))

ggsave(paste0(figDir, "spinup_init_matchedsites_bmdiff.png"))

mgmtdat = as.data.frame(mgmtdat)
mgmtdat$classification2 = ""
mgmtdat[mgmtdat$classification %in% c("BSP"),]$classification2 = "Black spruce"
mgmtdat[mgmtdat$classification %in% c("LPIN"),]$classification2 = "Pine"
mgmtdat[mgmtdat$classification %in% c("MIX"),]$classification2 = "Mixed forest"
mgmtdat[mgmtdat$classification %in% c("APIN","WSB","ASP","ABI","BIR","REP","OTH","WSPP","WSPS"),]$classification2 = "Other"

years = 1915:2100
#exclude sites initialized in 1983 from this
mgmt = mgmtdat[mgmtdat$year %in% years,]
mgmt = count(mgmt, classification2, iscen, cscen, year) %>% spread(classification2, n, fill = 0)
mgmt = reshape2::melt(mgmt, id.vars = c("iscen","cscen","year"))
mgmt$cscen = factor(mgmt$cscen, levels = c("hist","gcm45","gcm85"), labels = c("Historical climate","SSP2-4.5", "SSP5-8.5"))

# summarize initialization strength
temp = as.data.frame(dat %>% group_by(cscen, iscen, year) %>% summarize(signal = mean(init_strength, na.rm = T)))
temp$cscen = factor(temp$cscen, levels = c("hist","gcm45","gcm85"), labels = c("Historical climate","SSP2-4.5", "SSP5-8.5"))

mgmt = (merge(mgmt, temp, by = c("cscen","iscen","year")))
mgmt$iscen = factor(mgmt$iscen, levels = c("Spin-up","Initialization"))

g<-ggplot(mgmt) + facet_grid(iscen~cscen) + geom_line(aes(x = year, y = value, color = variable)) +  geom_line(aes(x = year, y = signal*500)) + 
  geom_vline(xintercept = c(2020,2060,2100), linetype = 3) + 
  scale_y_continuous("Initialization strength", sec.axis = sec_axis(~ . / 500, name = "Initialization strength")) + 
  ylab("Number of sites") + xlab("Year") + 
  scale_color_manual(breaks = c("Black spruce","Mixed forest","Pine","Other","Initialization strength"), 
                     values = c("#1F78BA","#B2DF8A","#33A02C","#FDBF6F","#000000"), name = "Forest type") +
  theme_linedraw() + theme(axis.text = element_text(size = 15))
ggsave(paste0(figDir, "spinup_init_classdiff.png"), g, height = 7, width = 14)

# Supplementary fig with 500-year spinup 
dat = fread("output_data/Landsat_init/spinup/500y/Active_Management_Data.csv")
dat$classification = sapply(dat$classification, FUN = function(x) classes11hash[[toString(x)]])
years = seq(0, 500, 5)
dat$classification = factor(dat$classification)
dat =  as.data.frame(dat %>% group_by(siteID) %>% subset(year %in% years) %>% mutate(nextclass = lead(classification, n = 1)))
dat = as.data.frame(dat %>% group_by(siteID) %>% mutate(nextyear = lead(year, n = 1)))

g <- ggplot(dat, aes(x = year, next_x = nextyear,
                     node = classification, next_node = nextclass,
                     fill = factor(classification))) +
  geom_sankey(flow.alpha = 0.9, width = 0.5) +
  theme_linedraw() + scale_fill_manual(breaks = classes11long, values = classpal11, name = "Forest classification") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(colour = 'black', size = 15),
        # legend.position = c(0.2,0.25),
        axis.text = element_text(size = 15)) +
  ylab("Number of sites") + xlab("Year") + ggtitle("UVAFME spinup under historic climate KLC")
ggsave(paste0(figDir, "spinup_500y_KLC.png"), g, width = 10, height = 6)

