#
# #########################
# Purpose: All figures and analyses for manuscript. Copying final versions for pub here so it's quick to update figs and stats if model is re-run
# Author: Shelby W. Sundquist
# Date: Feb, 2023
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
source("scripts/analyses/TVSF/Utilities.R")

figDir = "scripts/analyses/TVSF/figs/"

## First: check that all model runs have completed without errors. 
## Jobnames for any erroneous runs will be output in a txt and this script will be halted

### Checking landscape-level outputs
goodlastline = "   ============================================================================="
joberrors = c()
errorcodes = c()
for(m in amgmts[1]) {
  for(c in clims) {
    for(u in nmunits) {
      fname = paste0(logDir, landscapeDir, m, "_", c, "_", u, "_err.log")
      dat = readLines(fname)
      nsites = as.integer(str_extract(dat[6], "[0-9]+"))
      sitesrun = length( dat[grepl("Running for site", dat)])
      lastline = dat[length(dat)]
      if(lastline == goodlastline & nsites == sitesrun) {
        # scenario ran correctly
      } else {
        jobname = paste0("sbatch ",jobDir, "mgmt_testing/", m, "_", c, "_", u, ".slurm")
        joberrors = c(joberrors, jobname) 
        errorcodes = c(errorcodes, fname, lastline, paste0(sitesrun, " out of ", nsites, " run"))
      }
    }
  }
}

if(length(errorcodes > 0)) {
  badjobs = paste0(jobDir, "failedjobs.txt")
  fileConn <- file(badjobs)
  writeLines(joberrors, fileConn)
  close(fileConn)
  
  badjobs = paste0(jobDir, "failedjobserrors.txt")
  fileConn <- file(badjobs)
  writeLines(errorcodes, fileConn)
  close(fileConn)
  
  stop(paste0("Not all model runs have successfully completed. See failedjobs and failedjoberrors in ", jobDir))
}

## All model runs were successful: time to compile any new outputs
# This is very time consuming; don't run it every time 

fnames = c()
fnames = c("Total_Plot_Values.csv","Forestry_Data.csv","Active_Management_Data.csv","Climate.csv",
           "Dead_Species_Data.csv","SoilDecomp.csv","Fire_Summary.csv") #"Species_Data.csv",

# combine outputs for no management
for(f in fnames) {
  print(f)
  for(m in amgmts[1]) {
    for(c in clims) {
      dat = data.frame(matrix(ncol=length(cols[[f]]),nrow=0, dimnames=list(NULL, cols[[f]])))
      for(u in nmunits) { 
        #read all units for scenario 
        path = paste0(outDir, landscapeDir, m, "/", c, "/unit_", u, "/", f)
        read = fread(path)
        dat = rbind(dat, read)
      }
      # write file with all units combined for each scenario
      path = paste0(outDir, landscapeDir, m, "/", c, "/", f)
      fwrite(dat, path)
    }
  }
}

stop("Compiled files.")


## New outputs are now compiled: time to generate figs and stats
### 1. Initialization vs spinup forest types (sankey plot)
f = "Active_Management_Data.csv"
mgmtdat = data.frame(matrix(ncol=length(cols[[f]]) + 1,nrow=0, dimnames=list(NULL, c(cols[[f]], "iscen"))))


mgmttemp = fread(paste0(outDir, "mgmt_testing/cf/hist/unit_4C/",  f))
mgmttemp$iscen = "Initialization"
mgmtdat = rbind(mgmtdat, mgmttemp)

mgmttemp = fread(paste0(outDir, "inits_testing/spinup/", f))
mgmttemp$iscen = "Spin-up"
mgmtdat = rbind(mgmtdat, mgmttemp)

mgmtdat$classification = sapply(mgmtdat$classification, FUN = function(x) classes10hash[[toString(x)]])
mgmttemp = mgmtdat # snap of current dataset
years = seq(1913, 2100, 5)
mgmtdat$classification = factor(mgmtdat$classification)
mgmtdat$iscen = factor(mgmtdat$iscen, levels = c("Spin-up","Initialization"))
mgmtdat =  as.data.frame(mgmtdat %>% group_by(siteID, iscen) %>% subset(year %in% years) %>% mutate(nextclass = lead(classification, n = 1)))
mgmtdat = as.data.frame(mgmtdat %>% group_by(siteID,  iscen) %>% mutate(nextyear = lead(year, n = 1)))

g <- ggplot(subset(mgmtdat), aes(x = year, next_x = nextyear,
                                 node = classification, next_node = nextclass,
                                 fill = factor(classification))) +
  geom_sankey(flow.alpha = 0.9, width = 0.5) +
  facet_wrap(~iscen, nrow = 2) +
  theme_linedraw() + scale_fill_manual(breaks = classes10long, values = classpal10, name = "Forest classification") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(colour = 'black', size = 15),
        # legend.position = c(0.2,0.25),
        axis.text = element_text(size = 15)) +
  ylab("Number of sites") + xlab("Year") + ggtitle("UVAFME initialization vs spinup under historic climate")
ggsave(paste0(figDir, "init_spinup_hist_sankey_4C.png"), g, width = 10, height = 8)


g <- ggplot(subset(mgmtdat), aes(x = year, next_x = nextyear,
                                 node = classification, next_node = nextclass,
                                 fill = factor(classification))) +
  geom_sankey(flow.alpha = 0.9, width = 0.5) +
  facet_wrap(~iscen, nrow = 2) +
  scale_fill_manual(breaks = classes10long, values = classpal10, name = "Forest classification") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        strip.text = element_text(colour = 'black'),
        strip.background = element_rect(fill="transparent", color = NA)) +
  ylab("Number of sites") + xlab("Year") + ggtitle("UVAFME initialization vs spinup under historic climate") +
  transparent 
ggsave(paste0(figDir,"init_spinup_hist_sankey_trans.png"), g, width = 10, height = 8)

#### 1.2 Spinup vs init paired-sites comparison for unit 4C
f = "Active_Management_Data.csv"
mgmtdat = data.frame(matrix(ncol=length(cols[[f]]) + 2,nrow=0, dimnames=list(NULL, c(cols[[f]], "iscen", "cscen"))))


mgmttemp = fread(paste0(outDir, "inits_testing/spinup/", f))
mgmttemp$iscen = "Spin-up"; mgmttemp$cscen = "hist"
mgmtdat = rbind(mgmtdat, mgmttemp)
mgmttemp = fread(paste0(outDir, "inits_testing/spinup/r45/", f))
mgmttemp$iscen = "Spin-up"; mgmttemp$cscen = "gcm45"
mgmtdat = rbind(mgmtdat, mgmttemp)
mgmttemp = fread(paste0(outDir, "inits_testing/spinup/r85/", f))
mgmttemp$iscen = "Spin-up"; mgmttemp$cscen = "gcm85"
mgmtdat = rbind(mgmtdat, mgmttemp)

mgmttemp = fread(paste0(outDir, "mgmt_testing/cf/hist/unit_4C/",  f))
mgmttemp$iscen = "Initialization"; mgmttemp$cscen = "hist"
mgmtdat = rbind(mgmtdat, mgmttemp)
mgmttemp = fread(paste0(outDir, "mgmt_testing/cf/gcm45/unit_4C/",  f))
mgmttemp$iscen = "Initialization"; mgmttemp$cscen = "gcm45"
mgmtdat = rbind(mgmtdat, mgmttemp)
mgmttemp = fread(paste0(outDir, "mgmt_testing/cf/gcm85/unit_4C/",  f))
mgmttemp$iscen = "Initialization"; mgmttemp$cscen = "gcm85"
mgmtdat = rbind(mgmtdat, mgmttemp)

f = "Total_Plot_Values.csv"
dat = data.frame(matrix(ncol=length(cols[[f]]) + 2,nrow=0, dimnames=list(NULL, c(cols[[f]], "iscen", "cscen"))))

temp = fread(paste0(outDir, "inits_testing/spinup/", f))
temp$iscen = "Spin-up"; temp$cscen = "hist"
dat = rbind(dat, temp)
temp = fread(paste0(outDir, "inits_testing/spinup/r45/", f))
temp$iscen = "Spin-up"; temp$cscen = "gcm45"
dat = rbind(dat, temp)
temp = fread(paste0(outDir, "inits_testing/spinup/r85/", f))
temp$iscen = "Spin-up"; temp$cscen = "gcm85"
dat = rbind(dat, temp)
dat$init_strength = NA

temp = fread(paste0(outDir, "mgmt_testing/cf/hist/unit_4C/",  f))
temp$iscen = "Initialization"; temp$cscen = "hist"
dat = rbind(dat, temp)
temp = fread(paste0(outDir, "mgmt_testing/cf/gcm45/unit_4C/",  f))
temp$iscen = "Initialization"; temp$cscen = "gcm45"
dat = rbind(dat, temp)
temp = fread(paste0(outDir, "mgmt_testing/cf/gcm85/unit_4C/",  f))
temp$iscen = "Initialization"; temp$cscen = "gcm85"
dat = rbind(dat, temp)

dat = dat[order(dat$iscen, decreasing = F),]

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
mgmtdat[mgmtdat$classification %in% c("WSPP","WSPS"),]$classification2 = "White spruce"
mgmtdat[mgmtdat$classification %in% c("WSB","ASP","ABI","BIR","MIX","REP","OTH"),]$classification2 = "Deciduous, mixed, other"

mgmt = as.data.frame(mgmtdat %>% group_by(iscen, cscen, classification2) %>% subset(year == 2100) %>% reframe(n()))
years = 1913:2100
#exclude sites initialized in 1983 from this
mgmt = mgmtdat[mgmtdat$year %in% years & !mgmtdat$siteID %in% c(10868,10881,10891,10904,10920,10936),]
mgmt = count(mgmt, classification2, iscen, cscen, year) %>% spread(classification2, n, fill = 0)
mgmt = reshape2::melt(mgmt, id.vars = c("iscen","cscen","year"))
mgmt$cscen = factor(mgmt$cscen, levels = c("hist","gcm45","gcm85"), labels = c("Historical climate","SSP2-4.5", "SSP5-8.5"))

# summarize initialization strength
temp = as.data.frame(dat %>% group_by(cscen, iscen, year) %>% subset(!siteID %in% c(10868,10881,10891,10904,10920,10936)) %>% summarize(signal = mean(init_strength, na.rm = T)))
temp$cscen = factor(temp$cscen, levels = c("hist","gcm45","gcm85"), labels = c("Historical climate","SSP2-4.5", "SSP5-8.5"))

mgmt = (merge(mgmt, temp, by = c("cscen","iscen","year")))
mgmt$iscen = factor(mgmt$iscen, levels = c("Spin-up","Initialization"))

g<-ggplot(mgmt) + facet_grid(iscen~cscen) + geom_line(aes(x = year, y = value, color = variable)) +  geom_line(aes(x = year, y = signal*704)) + 
  geom_vline(xintercept = c(2020,2060,2100), linetype = 3) + 
  scale_y_continuous("Initialization strength", sec.axis = sec_axis(~ . / 704, name = "Initialization strength")) + 
  ylab("Number of sites") + xlab("Year") + 
  scale_color_manual(breaks = c("Black spruce","White spruce","Deciduous, mixed, other","Initialization strength"), values = c("#1F78B4","#FB9A99","#B2Df8A","#000000"), name = "Forest type") +
  theme_linedraw() + theme(axis.text = element_text(size = 15))
ggsave(paste0(figDir, "spinup_init_classdiff.png"), g, height = 7, width = 14)

## Remove large dataframes between tests to avoid confusion, save space
rm(dat)
rm(mgmttemp)
rm(mgmtdat)
rm(mgmt)


### 2. 20-year landscape-wide biomass density and forest classifications between climate scenarios

f = "Active_Management_Data.csv"
mgmtdat = data.frame(matrix(ncol=length(cols[[f]]) + 2,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen","mgmtscen"))))
for(m in amgmts) {
  for(c in clims) {
    path = paste0(outDir, landscapeDir, m, "/", c, "/", f)
    temp = fread(path)
    temp$mgmtscen = m
    temp$climscen = c
    mgmtdat = rbind(mgmtdat, temp)
  }
}

f = "Total_Plot_Values.csv"
plotsdat = data.frame(matrix(ncol=length(cols[[f]]) + 2,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen","mgmtscen"))))
for(m in amgmts) {
  for(c in clims) {
    path = paste0(outDir, landscapeDir, m, "/", c, "/", f)
    temp = fread(path)
    temp$mgmtscen = m
    temp$climscen = c
    plotsdat = rbind(plotsdat, temp)
  }
}

dat = merge(plotsdat, mgmtdat, by = c("siteID","year","mgmtscen","climscen"))
dat = as.data.frame(dat)
dat$classification2 = ""
dat[dat$classification == "BSP",]$classification2 = "Black spruce"
dat[dat$classification %in% c("WSPP","WSPS"),]$classification2 = "White spruce"
dat[dat$classification %in% c("ASP","ABI","BIR"),]$classification2 = "Deciduous (aspen/birch)"
dat[dat$classification %in% c("MIX","WSB"),]$classification2 = "Mixed forest"
dat[dat$classification %in% c("REP","OTH"),]$classification2 = "Reproductive/ other"

dat$climscen2 = climnames(as.character(dat$climscen))
dat$climscen2 = factor(dat$climscen2, c("Historic climate","Moderate climate warming","Extreme climate warming"))

years = c(2020,2040,2060,2080,2100)
g <- ggplot(subset(dat, year %in% years & mgmtscen == "cf")) + 
  geom_density(aes(x = total_biomC, fill = classification2, color = classification2,  alpha = 0.33), position = "stack") + 
  facet_grid(climscen2~year, switch = "y")  + ylim(0,0.1) + 
  theme_linedraw() + scale_fill_manual(breaks = classes5long, values = classpal5, name = "Forest type") +
  scale_color_manual (breaks = classes5long, values = classpal5, name = "Forest type") + 
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(colour = 'black')) +
  xlab("Biomass (tC/ha)") + ylab("Density") + scale_alpha(guide = 'none')

ggsave(paste0(figDir, "clim_biomass_class_20yrs.png"), g, width = 12, height = 8)


g <- ggplot(subset(dat, year %in% years & mgmtscen == "cf")) + 
  geom_density(aes(x = total_biomC, fill = classification2, color = classification2,  alpha = 0.33), position = "stack") + 
  facet_grid(climscen2~year, switch = "y")  + ylim(0,0.25) + 
  scale_fill_manual(breaks = classes5long, values = classpal5, name = "Forest type") +
  scale_color_manual (breaks = classes5long, values = classpal5, name = "Forest type") + 
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        strip.text = element_text(colour = 'black')) +
  xlab("Biomass (tC/ha)") + ylab("Density") + scale_alpha(guide = 'none') + transparent

ggsave(paste0(figDir, "clim_biomass_class_20yrs_trans.png"), g, width = 12, height = 8)


### 2.3.2 Baseline is historic climate run
years = c(2014, 2020,2040,2060,2080,2100)
#average biomass per site by forest type
foresttype_biomass_summary = matrix(nrow = 0, ncol = 6)
colnames(foresttype_biomass_summary) = c("classification", "mn_biomass","se_biomass","n", "year", "clim")

for(y in years) {
  for(c in clims) {
    ft_biomass = dat %>% group_by(classification2) %>% subset(year == y & mgmtscen == "cf" & climscen == c) %>% summarize(mean(total_biomC, na.rm = T))
    ft_biomass_se = dat %>% group_by(classification2) %>% subset(year == y & mgmtscen == "cf" & climscen == c) %>% summarize(sd(total_biomC, na.rm = T)/sqrt(n()))
    ft_sites = dat %>% group_by(classification2) %>% subset(year == y & mgmtscen == "cf" & climscen == c) %>% summarize(n())
    ft_biomass = cbind(ft_biomass, ft_biomass_se[,2])
    ft_biomass = cbind(ft_biomass, ft_sites[,2])
    colnames(ft_biomass) = c("classification", "mn_biomass","se_biomass","n")
    ft_biomass$clim = c
    ft_biomass$year = y
    foresttype_biomass_summary = rbind(foresttype_biomass_summary, ft_biomass)
  } 
}

#biomass difference
foresttype_biomass_summary$diff = NA
diff = foresttype_biomass_summary[foresttype_biomass_summary$clim == "gcm45",]$mn_biomass - foresttype_biomass_summary[foresttype_biomass_summary$clim == "hist",]$mn_biomass
foresttype_biomass_summary[foresttype_biomass_summary$clim == "gcm45",]$diff = diff
diff = foresttype_biomass_summary[foresttype_biomass_summary$clim == "gcm85",]$mn_biomass - foresttype_biomass_summary[foresttype_biomass_summary$clim == "hist",]$mn_biomass
foresttype_biomass_summary[foresttype_biomass_summary$clim == "gcm85",]$diff = diff

#relative biomass difference
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

fwrite(temp, paste0(figDir, "supp_t4.csv"))

g<- ggplot(subset(foresttype_biomass_summary, clim != "hist" & classification != "Reproductive/ other"), aes(x = year, y = rdiff, group = year,
                                                               ymax = rdiff+ 1.96*rse_biomass, ymin = rdiff - 1.96*rse_biomass)) +
  geom_pointrange(aes(color = classification, , fatten = n/800), position =position_dodge2(width = 7)) + geom_hline(yintercept = 1) +
  geom_point(shape = 21, aes(size = n/50), position = position_dodge2(width = 7)) +
  facet_wrap(~clim, nrow = 1) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(breaks = classes5long, values = classpal5, name = "Forest classification") +
  theme_linedraw() + ylab("Biomass change compared to historic climate") + xlab("Year") +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(colour = 'black'), 
        axis.text = element_text(size = 15)) + guides(size = "none")
ggsave(paste0(figDir, "exogenous_classes_biomass_clims.png"), g, width = 12, height = 4)

g<- ggplot(subset(foresttype_biomass_summary, clim != "hist" & classification != "Reproductive/ other"), aes(x = year, y = rdiff, group = year,
                                                               ymax = rdiff+ 1.96*rse_biomass, ymin = rdiff - 1.96*rse_biomass)) +
  geom_pointrange(aes(color = classification, , fatten = n/800), position =position_dodge2(width = 7)) + geom_hline(yintercept = 1) +
  geom_point(shape = 21, aes(size = n/50), position = position_dodge2(width = 7)) +
  facet_wrap(~clim, nrow = 1) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(breaks = classes5long, values = classpal5, name = "Forest classification") +
  theme_linedraw() + ylab("Biomass change compared to historic climate") + xlab("Year") +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(colour = 'black'), 
        axis.text = element_text(size = 15)) + guides(size = "none") + transparent
ggsave(paste0(figDir, "exogenous_classes_biomass_clims_transparent.png"), g, width = 12, height = 4)


#Landscape-wide biomass
biomass_summary = matrix(nrow = 0, ncol = 4)
colnames(biomass_summary) = c("biomass","se_biomass","year", "clim")

years = c(2020,2060,2100)
for(y in years) {
  temp = (dat %>% group_by(climscen) %>% subset(mgmtscen == "cf" & year == y) %>% summarize(sum(total_biomC*25))) #convert ha to site area (25ha)
  temp2 = (dat %>% group_by(climscen) %>% subset(mgmtscen == "cf" & year == y) %>% summarize(sqrt(mean((total_biomC_sd*25)^2))/sqrt(n())))
  temp = cbind(temp, temp2[,2])
  colnames(temp) = c("clim","biomass","se_biomass")
  temp$year = y 
  biomass_summary = rbind(biomass_summary, temp)
}
biomass_summary$clim = factor(biomass_summary$clim, levels = c("hist","gcm45","gcm85"), labels = c("Historical climate","SSP2-4.5", "SSP5-8.5"))

g<- ggplot(biomass_summary, aes(x = as.factor(year), y = biomass, fill = clim)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = biomass - 1.96*se_biomass, ymax = biomass + 1.96*se_biomass), position = position_dodge()) + 
  scale_fill_manual(breaks = c("Historical climate","SSP2-4.5", "SSP5-8.5"), values = c("#D1E5F0", "#F4A582", "#B2182B"), name = "Climate scenario") + 
  theme_linedraw() + theme(axis.text = element_text(size = 15)) + ylab("TVSF-wide Biomass (tC)") + xlab("Year")
  
ggsave(paste0(figDir, "totalbiomass_clims.png"), g)

g<- ggplot(biomass_summary, aes(x = as.factor(year), y = biomass, fill = clim)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = biomass - 1.96*se_biomass, ymax = biomass + 1.96*se_biomass), position = position_dodge()) + 
  scale_fill_manual(breaks = c("Historical climate","SSP2-4.5", "SSP5-8.5"), values = c("#D1E5F0", "#F4A582", "#B2182B"), name = "Climate scenario") + 
  theme_linedraw() + theme(axis.text = element_text(size = 15)) + ylab("TVSF-wide Biomass (tC)") + xlab("Year") + transparent
  
ggsave(paste0(figDir, "totalbiomass_clims_transparent.png"), g)

###2.4 

dat = mgmtdat[mgmtdat$mgmtscen == "cf",]
dat$classification2 = ""
dat[dat$classification == "BSP",]$classification2 = "Black spruce"
dat[dat$classification %in% c("WSPP","WSPS"),]$classification2 = "White spruce"
dat[dat$classification %in% c("ASP","ABI","BIR"),]$classification2 = "Deciduous (aspen/birch)"
dat[dat$classification %in% c("MIX","WSB"),]$classification2 = "Mixed forest"
dat[dat$classification %in% c("REP","OTH"),]$classification2 = "Reproductive/ other"

dat$classification = dat$classification2

years = c(2014, seq(2020, 2100, 5))
dat$classification = factor(dat$classification)
dat =  as.data.frame(dat %>% group_by(siteID, climscen) %>% subset(year %in% years) %>% mutate(nextclass = lead(classification, n = 1)))
dat = as.data.frame(dat %>% group_by(siteID, climscen) %>% mutate(nextyear = lead(year, n = 1)))


g <- ggplot(subset(dat), aes(x = year, next_x = nextyear,
                                 node = classification, next_node = nextclass,
                                 fill = factor(classification))) +
  facet_wrap(~climscen, nrow = 3) +
  geom_sankey(flow.alpha = 0.9, width = 0.5) +theme_linedraw() +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(colour = 'black'), axis.text = element_text(size = 15)) +  
  scale_fill_manual(breaks = classes5long, values = classpal5, name = "Forest classification") + 
  geom_vline(xintercept = c(2020,2040,2060,2080,2100), linetype = 3) +
  ylab("Number of sites") + xlab("Year") 
ggsave(paste0(figDir, "sankey_clims_broadclasses.png"), g, width = 10, height = 8)


g <- ggplot(subset(dat), aes(x = year, next_x = nextyear,
                                 node = classification, next_node = nextclass,
                                 fill = factor(classification))) +
  facet_wrap(~climscen, nrow = 3) +
  geom_sankey(flow.alpha = 0.9, width = 0.5) +theme_linedraw() +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(colour = 'black'), axis.text = element_text(size = 15)) +  
  scale_fill_manual(breaks = classes5long, values = classpal5, name = "Forest classification") + 
  geom_vline(xintercept = c(2020,2040,2060,2080,2100), linetype = 3) +
  ylab("Number of sites") + xlab("Year") + transparent 
ggsave(paste0(figDir, "sankey_clims_broadclasses_transparent.png"), g, width = 10, height = 8)

rm(mgmtdat)
rm(plotsdat)
rm(dat)

# Supplementary figures: drought affects deciduous productivity, warmth releases bsp from nutrient limitation 

f = "Species_Data.csv"
specdat = data.frame(matrix(ncol=length(cols[[f]]) + 1,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen"))))

m = amgmts[1]
for(c in clims) {
  path = paste0(outDir, landscapeDir, m, "/", c, "/", f)
  temp = fread(path)
  temp$climscen = c
  specdat = rbind(specdat, temp)
}

specdat[specdat == -999] <- NA

temp = reshape2::melt(specdat[, c("siteID","year","climscen","species","degday_resp","drought_resp","shade_resp","perm_resp","flood_resp","nutrient_resp")], id.vars = c("siteID","year","climscen","species"))
# temp = temp[temp$variable %in% c("degday_resp","drought_resp","shade_resp","perm_resp","flood_resp","nutrient_resp")]
temp = temp %>% group_by(climscen, year, variable, species) %>% summarize(resp = mean(value, na.rm = T), lower = quantile(value, 0.025, na.rm = T), upper = quantile(value, 0.975, na.rm = T))
temp = as.data.frame(temp)

temp$variable = factor(temp$variable, levels = c("degday_resp","drought_resp","shade_resp","perm_resp","flood_resp","nutrient_resp"), 
                            labels = c("Temperature","Drought","Shade","Permafrost","Flood","Nutrient"))

temp$climscen = factor(temp$climscen, levels = c("hist","gcm45","gcm85"), labels = c("Historical climate","SSP2-4.5", "SSP5-8.5"))


g <- ggplot(subset(temp, species == "POPUtrem" & year >= 2013), aes(x = year, y = resp, color = climscen, fill = climscen)) + facet_wrap(~variable) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) + geom_line() + guides(fill=guide_legend(title="Climate scenario"),color=guide_legend(title="Climate scenario")) +
  ylab("Response (0-1)") + xlab("Year") +  ggtitle("Growth-limiting factors for TVSF aspen")

ggsave(paste0(figDir, "Aspen_limiting_factors.png"),g,  width= 10, height = 4)



g <- ggplot(subset(temp, species == "PICEmari" & year >= 2013), aes(x = year, y = resp, color = climscen, fill = climscen)) + facet_wrap(~variable) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) + geom_line()  + guides(fill=guide_legend(title="Climate scenario"),color=guide_legend(title="Climate scenario")) +
  ylab("Response (0-1)") + xlab("Year") +  ggtitle("Growth-limiting factors for TVSF black spruce")

ggsave(paste0(figDir, "Bsp_limiting_factors.png"),g,  width= 10, height = 4)


# Supplementary figure for dwarf black spruce init conditions
dat = fread(paste0(outDir, "mgmt_testing/bsptest/Species_Data.csv"))
dat[dat == -999] <- NA

g <- ggplot(dat, aes(year, total_biomC, fill = species)) + geom_area() + facet_wrap(~siteID) + theme_linedraw() + 
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(colour = 'black')) +
  ylab("Biomass (tC/ha)") + xlab("Year")
ggsave(paste0(figDir, "bspinit.png"), g, width = 15, height = 12)


temp = reshape2::melt(dat[,c(1,3,5,6:15)], id.vars = c("siteID","year","species"))
ggplot(temp[temp$species == "PICEmari",], aes(x = year, y = value, color = variable)) + facet_wrap(~siteID) + geom_line()
+ geom_line(data = dat, aes(x = year, y = stem_coun))


# Supplementary fig with 500-year spinup 
dat = fread("output_data/TVSF/inits_testing/spinup/500y/Active_Management_Data.csv")
dat$classification = sapply(dat$classification, FUN = function(x) classes10hash[[toString(x)]])
years = seq(0, 500, 5)
dat$classification = factor(dat$classification)
dat =  as.data.frame(dat %>% group_by(siteID) %>% subset(year %in% years) %>% mutate(nextclass = lead(classification, n = 1)))
dat = as.data.frame(dat %>% group_by(siteID) %>% mutate(nextyear = lead(year, n = 1)))

g <- ggplot(dat, aes(x = year, next_x = nextyear,
                     node = classification, next_node = nextclass,
                     fill = factor(classification))) +
  geom_sankey(flow.alpha = 0.9, width = 0.5) +
  theme_linedraw() + scale_fill_manual(breaks = classes10long, values = classpal10, name = "Forest classification") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(colour = 'black', size = 15),
        # legend.position = c(0.2,0.25),
        axis.text = element_text(size = 15)) +
  ylab("Number of sites") + xlab("Year") + ggtitle("UVAFME spinup under historic climate TVSF")
ggsave(paste0(figDir, "spinup_500y_TVSF.png"), g, width = 10, height = 6)

