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
source("/home/ss3988/UVAFME/UVAFME_scripts/analyses/TVSF/Utilities.R")

figDir = "/home/ss3988/UVAFME/UVAFME_scripts/analyses/TVSF/figs/"

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

goodendyear = 2100
for(m in amgmts[2:3]) {
  for(c in clims) {
    for(u in munits) {
      fname = paste0(logDir,landscapeDir, m, "_", c, "_", u, "_err.log")
      dat = readLines(fname)
      endyear = dat[length(dat)-1]
      lastline = dat[length(dat)]
      if(lastline == goodlastline & str_contains(endyear, goodendyear)) {
        # scenario ran correctly
      } else {
        jobname = paste0("sbatch ", jobDir, "mgmt_testing/", m, "_", c, "_", u, ".slurm")
        joberrors = c(joberrors, jobname) 
        errorcodes = c(errorcodes, fname, lastline, paste0(endyear, " out of ", goodendyear, " run"))
      }      
    }
  }
}

# ### Checking site-level outputs
# for(f in forests[1]) {
#   for(m in cmgmts) {
#     for(c in lclims) {
#       fname = paste0(logDir, siteDir, f, "_", m, "_", c, "_err.log")
#       dat = readLines(fname)
#       nsites = as.integer(str_extract(dat[6], "[0-9]+"))
#       sitesrun = length( dat[grepl("Running for site", dat)])
#       lastline = dat[length(dat)]
#       if(lastline == goodlastline & nsites == sitesrun) {
#         # scenario ran correctly
#       } else {
#         jobname = paste0("sbatch ",jobDir, "site_mgmt/", m, "_", c, "_", u, ".slurm")
#         joberrors = c(joberrors, jobname) 
#         errorcodes = c(errorcodes, fname, lastline, paste0(sitesrun, " out of ", nsites, " run"))
#       }
#     }
#   }
# }
# 
# 
# for(f in forests[2:3]) {
#   for(m in cmgmts[1:3]) {
#     for(c in lclims) {
#       fname = paste0(logDir, siteDir, f, "_", m, "_", c, "_err.log")
#       dat = readLines(fname)
#       nsites = as.integer(str_extract(dat[6], "[0-9]+"))
#       sitesrun = length( dat[grepl("Running for site", dat)])
#       lastline = dat[length(dat)]
#       if(lastline == goodlastline & nsites == sitesrun) {
#         # scenario ran correctly
#       } else {
#         jobname = paste0("sbatch ",jobDir, "site_mgmt/", f, "_", m, "_", c, ".slurm")
#         joberrors = c(joberrors, jobname) 
#         errorcodes = c(errorcodes, fname, lastline, paste0(sitesrun, " out of ", nsites, " run"))
#       }
#     }
#   }
# }

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


## Second: Compile output files from landscape runs (site outputs will be called in their own; there are fewer and they're easier to iterate through)
# Well I finally used there/their/they're in a sentence together

# This is very time consuming; don't run it every time 

fnames = c()
fnames = c("Total_Plot_Values.csv","Forestry_Data.csv","Active_Management_Data.csv","Climate.csv",
           "Dead_Species_Data.csv","SoilDecomp.csv","Fire_Summary.csv") #"Species_Data.csv",

# combine outputs for no management
print("working on counterfactual")
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

# # combine outputs for active management; these incorporate cf scenarios for units where management doesn't occur for landscape-wide summary
# print("working on active management")
# for(f in fnames) {
#   for(m in amgmts[2:3]) {
#     for(c in clims) {
#       dat = data.frame(matrix(ncol=length(cols[[f]]),nrow=0, dimnames=list(NULL, cols[[f]])))
#       for(u in munits) {
#         path = paste0(outDir, landscapeDir, m, "/", c, "/unit_", u, "/", f)
#         read = fread(path)
#         dat = rbind(dat, read)
#       }
#       for(u in nmunits[!nmunits %in% munits]) {
#         path = paste0(outDir, landscapeDir, "cf/", c, "/unit_", u, "/", f)
#         read = fread(path)
#         dat = rbind(dat, read)
#       }
#       # write file with all units combined for each scenario
#       path = paste0(outDir, landscapeDir, m, "/", c, "/", f)
#       fwrite(dat, path)
#     }
#   }
# }

stop("Compiled files.")



## New outputs are now compiled: time to generate figs and stats
### 1. Initialization vs spinup forest types (sankey plot)
f = "Active_Management_Data.csv"
mgmtdat = data.frame(matrix(ncol=length(cols[[f]]) + 1,nrow=0, dimnames=list(NULL, c(cols[[f]], "iscen"))))

# mgmttemp = fread(paste0(outDir, landscapeDir,  "cf/hist/unit_6/", f))
# mgmttemp = fread(paste0(outDir, "inits_testing/init/",  f))
mgmttemp = fread(paste0(outDir, "mgmt_testing/cf/hist/unit_4C/",  f))
mgmttemp$iscen = "Initialization"
mgmtdat = rbind(mgmtdat, mgmttemp)
# mgmttemp = fread(paste0(outDir, landscapeDir, "cf/hist/unit_6/spinup/", f))
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
ggsave("/home/ss3988/UVAFME/UVAFME_scripts/analyses/TVSF/figs/init_spinup_hist_sankey_4C.png", g, width = 10, height = 8)


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
ggsave("/home/ss3988/UVAFME/UVAFME_scripts/analyses/TVSF/figs/init_spinup_hist_sankey_trans.png", g, width = 10, height = 8)

#### 1.1 Forest classifications on a single-site species-level biomass vs time plot

#### Plot of wsp --> bsp
# mgmtdat = mgmttemp
# wsps = mgmtdat[mgmtdat$classification == "White spruce saw",]$siteID
# wsp_bsp = mgmtdat %>% group_by(siteID, climscen, iscen) %>% filter(siteID %in% wsps & classification == "Black spruce")
# specdat = fread(paste0(outDir, initDir, wsp_bsp$iscen[1],"/", wsp_bsp$climscen[1], "Species_Data.csv"))
# specdat = specdat[specdat$siteID == wsp_bsp$siteID[1],]
# specdat[specdat == -999] <- NA
# mgmt = subset(mgmtdat, siteID == wsp_bsp$siteID[1] & climscen == wsp_bsp$climscen[1] & iscen == wsp_bsp$iscen[1])
# mgmt = mgmt %>% group_by(classification) %>% summarize(start = min(year), end = max(year))
# 
# ggplot(specdat, aes(year, total_biomC, fill = species)) + geom_area() +
#   geom_vline(mgmt, xintercept = end) + geom_text(mgmt, aes(label = classification, x = (start + end)/2, y = 60)) + 
#   theme_linedraw() + scale_fill_manual(breaks = specs, values = specpal, name = "Species")

# I have nice versions of these from my comps prep!


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


# ggplot(endo_bm_diff, aes(x = cscen, y = diff)) + geom_boxplot() + geom_hline(yintercept = 0)

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

# ggplot(temp, aes(x = year, y = signal)) + facet_grid(iscen~cscen) + geom_line()

mgmt = (merge(mgmt, temp, by = c("cscen","iscen","year")))
mgmt$iscen = factor(mgmt$iscen, levels = c("Spin-up","Initialization"))

# mgmt[mgmt$iscen == "Initialization" & mgmt$year < 2013,]$value = 

g<-ggplot(mgmt) + facet_grid(iscen~cscen) + geom_line(aes(x = year, y = value, color = variable)) +  geom_line(aes(x = year, y = signal*704)) + 
  geom_vline(xintercept = c(2020,2060,2100), linetype = 3) + 
  scale_y_continuous("Initialization strength", sec.axis = sec_axis(~ . / 704, name = "Initialization strength")) + 
  ylab("Number of sites") + xlab("Year") + 
  scale_color_manual(breaks = c("Black spruce","White spruce","Deciduous, mixed, other","Initialization strength"), values = c("#1F78B4","#FB9A99","#B2Df8A","#000000"), name = "Forest type") +
  theme_linedraw() + theme(axis.text = element_text(size = 15))
ggsave(paste0(figDir, "spinup_init_classdiff.png"), g, height = 7, width = 14)
 
# mgmt = mgmt[order(mgmt$cscen, mgmt$classification2),]
# mgmt$`n()` = as.numeric(mgmt$`n()`)
# endo_class_diff = mgmt[mgmt$iscen == "Initialization",]$`n()` - mgmt[mgmt$iscen == "Spin-up",]$`n()`
# endo_class_diff = cbind(mgmt[mgmt$iscen == "Spin-up",], endo_class_diff)
# #how much more/less of each forest type does spinup have compared to init
# ggplot(endo_class_diff, aes(x = classification2, y = endo_class_diff, fill = classification2)) + 
#   geom_col(position = "dodge") + facet_wrap(~cscen, nrow = 1) + 
#   theme_linedraw() + theme(strip.background = element_rect(fill="white"), strip.text = element_text(colour = 'black')) + 
#   geom_hline(yintercept = 0)
#   # ylab("Forest type change compared to historic climate") + xlab("Year")
# 
# ggplot(mgmt, aes(y = `n()`, x = classification2, fill = iscen)) + facet_wrap(~cscen) + geom_col(position = "dodge")

## Remove large dataframes between tests to avoid confusion, save space
rm(dat)
rm(mgmttemp)
# rm(specdat)

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

# mgmtdat$mgmtscen = factor(mgmtdat$mgmtscen, amgmts)
# mgmtdat$climscen = factor(mgmtdat$climscen, clims)
# mgmtdat$year = as.integer(mgmtdat$year)

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
# plotsdat$climscen = factor(plotsdat$climscen, clims)
# plotsdat$mgmtscen = factor(plotsdat$mgmtscen, amgmts)


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


# 
# #### 2.1 Chi-sq statistics on landscape-wide management scenarios and forest classifications
# ##### Might add this later
# years = c(2060,2100)
# 
# #for forest classification
# for(y in years) {
#   cf_hist = summary(factor(dat[dat$mgmtscen == "cf" & dat$climscen == "hist" & dat$year == y,]$classification, levels = classes10))
#   cf_gcm45 = summary(factor(dat[dat$mgmtscen == "cf" & dat$climscen == "gcm45" & dat$year == y,]$classification, levels = classes10))
#   cf_gcm85 = summary(factor(dat[dat$mgmtscen == "cf" & dat$climscen == "gcm85" & dat$year == y,]$classification, levels = classes10))
#   bau_hist = summary(factor(dat[dat$mgmtscen == "bau" & dat$climscen == "hist" & dat$year == y,]$classification, levels = classes10))
#   bau_gcm45 = summary(factor(dat[dat$mgmtscen == "bau" & dat$climscen == "gcm45" & dat$year == y,]$classification, levels = classes10))
#   bau_gcm85 = summary(factor(dat[dat$mgmtscen == "bau" & dat$climscen == "gcm85" & dat$year == y,]$classification, levels = classes10))
#   prj_hist = summary(factor(dat[dat$mgmtscen == "prj" & dat$climscen == "hist" & dat$year == y,]$classification, levels = classes10))
#   prj_gcm45 = summary(factor(dat[dat$mgmtscen == "prj" & dat$climscen == "gcm45" & dat$year == y,]$classification, levels = classes10))
#   prj_gcm85 = summary(factor(dat[dat$mgmtscen == "prj" & dat$climscen == "gcm85" & dat$year == y,]$classification, levels = classes10))
#   
#   r1c1 = chisq.test(bau_hist, y = cf_hist)$p.value
#   r2c1 = chisq.test(prj_hist, y = cf_hist)$p.value
#   r3c1 = chisq.test(prj_hist, y = bau_hist)$p.value
#   
#   r1c2 = chisq.test(bau_gcm45, y = cf_gcm45)$p.value
#   r2c2 = chisq.test(prj_gcm45, y = cf_gcm45)$p.value
#   r3c2 = chisq.test(prj_gcm45, y = bau_gcm45)$p.value
#   
#   r1c3 = chisq.test(bau_gcm85, y = cf_gcm85)$p.value
#   r2c3 = chisq.test(prj_gcm85, y = cf_gcm85)$p.value
#   r3c3 = chisq.test(prj_gcm85, y = bau_gcm85)$p.value
#   
#   r1 = c(r1c1, r1c2, r1c3)
#   r2 = c(r2c1, r2c2, r2c3)
#   r3 = c(r3c1, r3c2, r3c3)
#   chisqtable = rbind(r1,r2, r3)
#   colnames(chisqtable) = c("Historical climate","SSP2-4.5","SSP5-8.5")
#   rownames(chisqtable) = c("BAU vs CF","ADP vs CF","BAU vs ADP")
#   chisqtable = round(chisqtable, 3)
#   write.csv(chisqtable, paste0(figDir,"2/forestclasses_chisq_table_",y,".csv"))
# }
# 
# #for biomass
# for(y in years) {
#   cf_hist = (dat[dat$mgmtscen == "cf" & dat$climscen == "hist" & dat$year == y,]$total_biomC)
#   cf_gcm45 = (dat[dat$mgmtscen == "cf" & dat$climscen == "gcm45" & dat$year == y,]$total_biomC)
#   cf_gcm85 = (dat[dat$mgmtscen == "cf" & dat$climscen == "gcm85" & dat$year == y,]$total_biomC)
#   bau_hist = (dat[dat$mgmtscen == "bau" & dat$climscen == "hist" & dat$year == y,]$total_biomC)
#   bau_gcm45 = (dat[dat$mgmtscen == "bau" & dat$climscen == "gcm45" & dat$year == y,]$total_biomC)
#   bau_gcm85 = (dat[dat$mgmtscen == "bau" & dat$climscen == "gcm85" & dat$year == y,]$total_biomC)
#   prj_hist = (dat[dat$mgmtscen == "prj" & dat$climscen == "hist" & dat$year == y,]$total_biomC)
#   prj_gcm45 = (dat[dat$mgmtscen == "prj" & dat$climscen == "gcm45" & dat$year == y,]$total_biomC)
#   prj_gcm85 = (dat[dat$mgmtscen == "prj" & dat$climscen == "gcm85" & dat$year == y,]$total_biomC)
#   
#   r1c1 = chisq.test(bau_hist, y = cf_hist)$p.value
#   r2c1 = chisq.test(prj_hist, y = cf_hist)$p.value
#   r3c1 = chisq.test(prj_hist, y = bau_hist)$p.value
#   
#   r1c2 = chisq.test(bau_gcm45, y = cf_gcm45)$p.value
#   r2c2 = chisq.test(prj_gcm45, y = cf_gcm45)$p.value
#   r3c2 = chisq.test(prj_gcm45, y = bau_gcm45)$p.value
#   
#   r1c3 = chisq.test(bau_gcm85, y = cf_gcm85)$p.value
#   r2c3 = chisq.test(prj_gcm85, y = cf_gcm85)$p.value
#   r3c3 = chisq.test(prj_gcm85, y = bau_gcm85)$p.value
#   
#   r1 = c(r1c1, r1c2, r1c3)
#   r2 = c(r2c1, r2c2, r2c3)
#   r3 = c(r3c1, r3c2, r3c3)
#   chisqtable = rbind(r1,r2, r3)
#   colnames(chisqtable) = c("Historical climate","SSP2-4.5","SSP5-8.5")
#   rownames(chisqtable) = c("BAU vs CF","ADP vs CF","BAU vs ADP")
#   chisqtable = round(chisqtable, 3)
#   write.csv(chisqtable, paste0(figDir,"2/biomass_chisq_table_",y,".csv"))
# }
# 
# #for fire
# f = "Fire_Summary.csv"
# firedat = data.frame(matrix(ncol=length(cols[[f]]) + 2,nrow=0, dimnames=list(NULL, c(cols[[f]], "mgmtscen","climscen"))))
# 
# for(m in amgmts) {
#   for(c in clims) {
#     path = paste0(outDir, landscapeDir, m, "/", c, "/", f)
#     temp = fread(path)
#     temp$mgmtscen = m
#     temp$climscen = c
#     firedat = rbind(firedat, temp)
#   }
# }
# dat = merge(dat, firedat, by = c("siteID","year","mgmtscen","climscen"))
# 
# #for fire probability
# for(y in years) {
#   cf_hist = (dat[dat$mgmtscen == "cf" & dat$climscen == "hist" & dat$year == y,]$NFires)
#   cf_gcm45 = (dat[dat$mgmtscen == "cf" & dat$climscen == "gcm45" & dat$year == y,]$NFires)
#   cf_gcm85 = (dat[dat$mgmtscen == "cf" & dat$climscen == "gcm85" & dat$year == y,]$NFires)
#   bau_hist = (dat[dat$mgmtscen == "bau" & dat$climscen == "hist" & dat$year == y,]$NFires)
#   bau_gcm45 = (dat[dat$mgmtscen == "bau" & dat$climscen == "gcm45" & dat$year == y,]$NFires)
#   bau_gcm85 = (dat[dat$mgmtscen == "bau" & dat$climscen == "gcm85" & dat$year == y,]$NFires)
#   prj_hist = (dat[dat$mgmtscen == "prj" & dat$climscen == "hist" & dat$year == y,]$NFires)
#   prj_gcm45 = (dat[dat$mgmtscen == "prj" & dat$climscen == "gcm45" & dat$year == y,]$NFires)
#   prj_gcm85 = (dat[dat$mgmtscen == "prj" & dat$climscen == "gcm85" & dat$year == y,]$NFires)
#   
#   r1c1 = chisq.test(bau_hist, y = cf_hist)$p.value
#   r2c1 = chisq.test(prj_hist, y = cf_hist)$p.value
#   r3c1 = chisq.test(prj_hist, y = bau_hist)$p.value
#   
#   r1c2 = chisq.test(bau_gcm45, y = cf_gcm45)$p.value
#   r2c2 = chisq.test(prj_gcm45, y = cf_gcm45)$p.value
#   r3c2 = chisq.test(prj_gcm45, y = bau_gcm45)$p.value
#   
#   r1c3 = chisq.test(bau_gcm85, y = cf_gcm85)$p.value
#   r2c3 = chisq.test(prj_gcm85, y = cf_gcm85)$p.value
#   r3c3 = chisq.test(prj_gcm85, y = bau_gcm85)$p.value
#   
#   r1 = c(r1c1, r1c2, r1c3)
#   r2 = c(r2c1, r2c2, r2c3)
#   r3 = c(r3c1, r3c2, r3c3)
#   chisqtable = rbind(r1,r2, r3)
#   colnames(chisqtable) = c("Historical climate","SSP2-4.5","SSP5-8.5")
#   rownames(chisqtable) = c("BAU vs CF","ADP vs CF","BAU vs ADP")
#   chisqtable = round(chisqtable, 3)
#   write.csv(chisqtable, paste0(figDir,"2/fireprob_chisq_table_",y,".csv"))
# }
# 
# 
# #for fire intensity
# for(y in years) {
#   cf_hist = (dat[dat$mgmtscen == "cf" & dat$climscen == "hist" & dat$year == y,]$I_r)
#   cf_gcm45 = (dat[dat$mgmtscen == "cf" & dat$climscen == "gcm45" & dat$year == y,]$I_r)
#   cf_gcm85 = (dat[dat$mgmtscen == "cf" & dat$climscen == "gcm85" & dat$year == y,]$I_r)
#   bau_hist = (dat[dat$mgmtscen == "bau" & dat$climscen == "hist" & dat$year == y,]$I_r)
#   bau_gcm45 = (dat[dat$mgmtscen == "bau" & dat$climscen == "gcm45" & dat$year == y,]$I_r)
#   bau_gcm85 = (dat[dat$mgmtscen == "bau" & dat$climscen == "gcm85" & dat$year == y,]$I_r)
#   prj_hist = (dat[dat$mgmtscen == "prj" & dat$climscen == "hist" & dat$year == y,]$I_r)
#   prj_gcm45 = (dat[dat$mgmtscen == "prj" & dat$climscen == "gcm45" & dat$year == y,]$I_r)
#   prj_gcm85 = (dat[dat$mgmtscen == "prj" & dat$climscen == "gcm85" & dat$year == y,]$I_r)
#   
#   r1c1 = chisq.test(bau_hist, y = cf_hist)$p.value
#   r2c1 = chisq.test(prj_hist, y = cf_hist)$p.value
#   r3c1 = chisq.test(prj_hist, y = bau_hist)$p.value
#   
#   r1c2 = chisq.test(bau_gcm45, y = cf_gcm45)$p.value
#   r2c2 = chisq.test(prj_gcm45, y = cf_gcm45)$p.value
#   r3c2 = chisq.test(prj_gcm45, y = bau_gcm45)$p.value
#   
#   r1c3 = chisq.test(bau_gcm85, y = cf_gcm85)$p.value
#   r2c3 = chisq.test(prj_gcm85, y = cf_gcm85)$p.value
#   r3c3 = chisq.test(prj_gcm85, y = bau_gcm85)$p.value
#   
#   r1 = c(r1c1, r1c2, r1c3)
#   r2 = c(r2c1, r2c2, r2c3)
#   r3 = c(r3c1, r3c2, r3c3)
#   chisqtable = rbind(r1,r2, r3)
#   colnames(chisqtable) = c("Historical climate","SSP2-4.5","SSP5-8.5")
#   rownames(chisqtable) = c("BAU vs CF","ADP vs CF","BAU vs ADP")
#   chisqtable = round(chisqtable, 3)
#   write.csv(chisqtable, paste0(figDir,"2/fireintensity_chisq_table_",y,".csv"))
# }
# 
# 

# 
# #### 2.2 20-year active management landscape biomass density plots; compare management scenarios
# years = c(2020,2040,2060,2080,2100)
# 
# plotsdat$mgmtscen2 = mgmtnames(as.character(plotsdat$mgmtscen))
# plotsdat$climscen2 = climnames(as.character(plotsdat$climscen))
# plotsdat$climscen2 = factor(plotsdat$climscen2, c("Historic climate","Moderate climate warming","Extreme climate warming"))
# plotsdat$mgmtscen2 = factor(plotsdat$mgmtscen2, c("No management","Business as usual","Projected"))
# 
# g <- ggplot(subset(plotsdat, year %in% years), aes(x = total_biomC, fill = mgmtscen2, alpha = 0.33)) + 
#   geom_density() + facet_grid(climscen2~year, switch = "y")  + 
#   theme_linedraw() + 
#   theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
#         strip.background = element_rect(fill="white"),
#         strip.text = element_text(colour = 'black')) +
#   xlab("Biomass (tC/ha)") + ylab("Density") + scale_alpha(guide = 'none') + guides(fill=guide_legend(title="Management scenario"), color = "none") 
# 
# ggsave(paste0(figDir, "mgmt_clim_biomass_density.png"), g, width = 12, height = 8)
# 


### 2.3 20-year relative changes in forest composition, biomass by forest type, and landscape-wide biomass
years = c(2014,2020,2040,2060,2080,2100)

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

foresttype_biomass_summary = foresttype_biomass_summary[order(foresttype_biomass_summary$year, foresttype_biomass_summary$clim, foresttype_biomass_summary$classification),]

base = foresttype_biomass_summary[foresttype_biomass_summary$year == 2014 & foresttype_biomass_summary$clim == "gcm45",]
foresttype_biomass_summary = foresttype_biomass_summary[foresttype_biomass_summary$year > 2014,]
foresttype_biomass_summary$rdiff = foresttype_biomass_summary$mn_biomass / rep(base$mn_biomass, 15)
foresttype_biomass_summary$rse_biomass = foresttype_biomass_summary$se_biomass / rep(base$mn_biomass, 15)
foresttype_biomass_summary$ndiff = foresttype_biomass_summary$n - rep(base$n, 15)

foresttype_biomass_summary$classification = factor(foresttype_biomass_summary$classification)

g<- ggplot(subset(foresttype_biomass_summary, clim != "hist" & classification != "Reproductive/ other"), aes(x = year, y = rdiff, group = year,
                                                               ymax = rdiff+ 1.96*rse_biomass, ymin = rdiff - 1.96*rse_biomass)) +
  geom_pointrange(aes(color = classification, fatten = n/800), position =position_dodge2(width = 7)) + geom_hline(yintercept = 1) +
  geom_point(shape = 21, aes(size = n/50), position = position_dodge2(width = 7)) +
  facet_wrap(~clim, nrow = 1) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(breaks = classes5long[1:4], values = classpal5[1:4], name = "Forest classification") +
  theme_linedraw() + ylab("Biomass change compared to historic climate") + xlab("Year") +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(colour = 'black')) + guides(size = "none")
ggsave(paste0(figDir, "2014baseline_classes_biomass_frequency_clims.png"), g, width = 12, height = 4)



# ggplot(subset(foresttype_biomass_summary, clim != "hist"), aes(x = year, y = rdiff, group = year,color = classification,
#                                                                ymax = rdiff+ 1.96*rse_biomass, ymin = rdiff - 1.96*rse_biomass)) + 
#   geom_pointrange(position = position_dodge2(width = 7) ) + geom_hline(yintercept = 1) +
#   facet_wrap(~clim, nrow = 1) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
#   scale_color_manual(breaks = classes5long, values = classpal5, name = "Forest classification") + 
#   theme_linedraw() + ylab("Biomass change compared to historic climate") + xlab("Year") + 
#   theme(strip.background = element_rect(fill="white"),
#         strip.text = element_text(colour = 'black'))

# 
# ggplot(subset(foresttype_biomass_summary, clim != "hist"), aes(x = as.factor(year), y = ndiff, fill = classification)) + 
#   geom_col(position = "dodge") + 
#   scale_fill_manual(breaks = classes5long, values = classpal5, name = "Forest classification") + facet_wrap(~clim, nrow = 1) + 
#   theme_linedraw() + theme(strip.background = element_rect(fill="white"), strip.text = element_text(colour = 'black')) + 
#   ylab("Forest type change compared to historic climate") + xlab("Year")
# ggsave(paste0(figDir, "2/classes_count_clims.png"), width = 12, height = 4)

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

# ggplot(subset(foresttype_biomass_summary, clim != "hist"), aes(x = year, y = rdiff, group = year,color = classification,
#                                                                ymax = rdiff+ 1.96*rse_biomass, ymin = rdiff - 1.96*rse_biomass)) + 
#   geom_pointrange(position = position_dodge2(width = 7) ) + geom_hline(yintercept = 1) +
#   facet_wrap(~clim, nrow = 1) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
#   scale_color_manual(breaks = classes5long, values = classpal5, name = "Forest classification") + 
#   theme_linedraw() + ylab("Biomass change compared to historic climate") + xlab("Year") + 
#   theme(strip.background = element_rect(fill="white"),
#         strip.text = element_text(colour = 'black'))


ggplot(subset(foresttype_biomass_summary, clim != "hist"), aes(x = as.factor(year), y = ndiff/15819, fill = classification)) + 
  geom_col(position = "dodge") + 
  scale_fill_manual(breaks = classes5long, values = classpal5, name = "Forest classification") + facet_wrap(~clim, nrow = 1) + 
  theme_linedraw() + theme(strip.background = element_rect(fill="white"), strip.text = element_text(colour = 'black'), 
                           axis.text = element_text(size = 15)) + 
  ylab("Forest type change compared to historic climate") + xlab("Year")
ggsave(paste0(figDir, "classes_count_clims.png"), width = 12, height = 4)


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
rm(firedat)
rm(dat)

### 3. Site-level effects on fire probability and intensity 

fire = data.frame(matrix(ncol=10,nrow=0, 
                         dimnames=list(NULL, c(cols[["Fire_Summary.csv"]], "mgmt","forest","clim"))))

clim = data.frame(matrix(ncol=23,nrow=0,
                         dimnames=list(NULL, c(cols[["Climate.csv"]],"mgmt","forest","clim"))))


mgmt = data.frame(matrix(ncol=7,nrow=0, 
                         dimnames=list(NULL, c(cols[["Active_Management_Data.csv"]],"mgmt","forest","clim"))))

cmgmts =c("cf","harv","thin30","shear","prune")
f = forests[1]
for(m in cmgmts) {
  for(c in lclims) {
    dir = paste0(f, "/", m,"/",c)
    temp = fread(paste0(outDir, siteDir, dir, "/Fire_Summary.csv"))
    temp$mgmt = m
    temp$forest = f
    temp$clim = c
    fire = rbind(fire, temp)
    
    temp = fread(paste0(outDir, siteDir, dir, "/Climate.csv"))
    temp$mgmt = m
    temp$forest = f
    temp$clim = c
    clim = rbind(clim, temp)
    
    temp = fread(paste0(outDir, siteDir, dir, "/Active_Management_Data.csv"))
    temp$mgmt = m
    temp$forest = f
    temp$clim = c
    mgmt = rbind(mgmt, temp)
  }
}

for(f in forests[2:3]) {
  for(m in cmgmts[1:3]) {
    for(c in lclims) {
      dir = paste0(f, "/", m,"/",c)
      temp = fread(paste0(outDir, siteDir, dir, "/Fire_Summary.csv"))
      temp$mgmt = m
      temp$forest = f
      temp$clim = c
      fire = rbind(fire, temp)
      
      temp = fread(paste0(outDir, siteDir, dir, "/Climate.csv"))
      temp$mgmt = m
      temp$forest = f
      temp$clim = c
      clim = rbind(clim, temp)
      
      temp = fread(paste0(outDir, siteDir, dir, "/Active_Management_Data.csv"))
      temp$mgmt = m
      temp$forest = f
      temp$clim = c
      mgmt = rbind(mgmt, temp)    
    }
  }
}

sites = fread(paste0(inDir, siteDir, cmgmts[1], "/UVAFME2018_site.csv"))
sites$mgmt = cmgmts[1]
for(m in cmgmts[2:5]) {
  temp = fread(paste0(inDir, siteDir, m, "/UVAFME2018_site.csv"))
  temp$mgmt = m
  sites = rbind(sites, temp) 
}

sites = assignforest(sites)
sites = sites[,c("site","mgmt","forest","year_management")]
sites$year_management = floor(sites$year_management)

fire = merge(fire, sites, by.x = c("siteID","forest","mgmt"), by.y = c("site","forest","mgmt"))
fire[fire$year_management == 0,]$year_management = 2102 #need this to be scheduled in future
fire$ysincemgmt = fire$year - fire$year_management
fire$ysincemgmt = sapply(fire$ysincemgmt, function(x) max(x, -2))


#### 3.1 Treatment ~ forest type time series
dat = merge(fire, mgmt, by = c("siteID","year","mgmt","forest","clim"))
dat = merge(dat, clim, by = c("siteID","year","mgmt","forest","clim"))

dat$classification2 = ""
dat[dat$classification == "BSP",]$classification2 = "Black spruce"
dat[dat$classification %in% c("WSPP","WSPS"),]$classification2 = "White spruce"
dat[dat$classification %in% c("ASP","ABI","BIR"),]$classification2 = "Deciduous (aspen/birch)"
dat[dat$classification %in% c("MIX","WSB"),]$classification2 = "Mixed forest"
dat[dat$classification %in% c("REP","OTH"),]$classification2 = "Reproductive/ other"


tdat = dat[dat$mgmt != "cf",]
tdat$classification = tdat$classification2
# tdat[tdat$mgmt == "cf",]$ysincemgmt = tdat[tdat$mgmt == "cf",]$year - 2014
# tdat$classification = sapply(tdat$classification, FUN = function(x) classes10hash[[toString(x)]])
tdat$classification = factor(tdat$classification)

tdat = tdat %>% group_by(siteID, mgmt, clim) %>% mutate(mndrydays = lead(rollmean(drydays, 5, fill = NA),2))
# tdat$drydaysgroup = (floor(tdat$mndrydays*20)+1)/20
tdat$drydaysgroup = round(tdat$mndrydays, digits = 1)
tdat = tdat %>% group_by(siteID, mgmt, clim) %>% mutate(drydaysgroup = if_else(ysincemgmt == 0, drydaysgroup, NA))
tdat = tdat %>% group_by(siteID, mgmt, clim) %>% mutate(drydaysgroup = max(drydaysgroup, na.rm = T))
# tdat = as.data.frame(tdat %>% group_by(siteID, mgmt, clim) %>% mutate(mgmtdrydays = max(mgmtdrydays, na.rm = T)))
# # tdat$mgmtdrydays = paste0(tdat$mgmtdrydays - 0.2, " - ", tdat$mgmtdrydays, " dry days")
# tdat[is.na(tdat$drydaysgroup),]$drydaysgroup = -1
# tdat = as.data.frame(tdat)
# tdat$mgmtdrydays = as.character(tdat$mgmtdrydays)
# tdat[!tdat$mgmtdrydays %in% c("0","0.05","0.1, 0.15, 0.2"),]$mgmtdrydays = "0.2-0.6"
# tdat[tdat$mgmtdrydays %in% c("0","0.05"),]$mgmtdrydays = "0-0.05"
# tdat[tdat$mgmtdrydays %in% c("0.1", "0.15", "0.2")]$mgmtdrydays = "0.05-0.0.2"
tdat$mgmtdrydays = ""
tdat[tdat$drydaysgroup == 0,]$mgmtdrydays = "< 0.1"
tdat[tdat$drydaysgroup == 0.1,]$mgmtdrydays = "0.1 - 0.2"
tdat[tdat$drydaysgroup >= 0.2,]$mgmtdrydays = "> 0.2"
tdat$mgmtdrydays = factor(tdat$mgmtdrydays, levels = (c("< 0.1","0.1 - 0.2","> 0.2")))

tdat = subset(tdat, ysincemgmt %in% c(1,5,10,20,30,40)) 
# year == -1 might not be proper forest type for management; management is prescribed for these sites
tdat = as.data.frame(tdat %>% group_by(siteID, mgmt, clim) %>% mutate(nextclass = lead(classification, n = 1)))
tdat = as.data.frame(tdat %>% group_by(siteID, mgmt, clim) %>% mutate(nextysince = lead(ysincemgmt, n = 1)))
tdat$mgmt = treatmentnames(tdat$mgmt)
tdat = tdat %>% group_by(ysincemgmt, mgmt, mgmtdrydays) %>% mutate(size = n())
tdat = tdat[order(abs(tdat$ysincemgmt)),]
tdat = tdat %>% group_by(mgmt, mgmtdrydays) %>% mutate(size = ifelse(ysincemgmt == 0, size, first(size)))

g <- ggplot(subset(tdat, mgmtdrydays == "< 0.1"), aes(x = ysincemgmt, next_x = nextysince, node = classification, next_node = nextclass, fill = factor(classification))) +
  geom_sankey(flow.alpha = 0.9, width = 0.5) + 
  facet_wrap(~mgmt, nrow = 1) + theme_linedraw() + scale_fill_manual(breaks = classes5long, values = classpal5, name = "Forest classification") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(colour = 'black')) + xlab ("Years since management") + ylab("Site-level forest classification") + 
  ggtitle("Climate impacts on forest transitions after management") + 
  geom_text(aes(label =  paste0( "n = ", size), x = Inf, y = -Inf, hjust=1, vjust=0))

ggsave("/home/ss3988/UVAFME/mgmt_results/figs/3/sitemgmt_drydays_least_sankey.png", g, width = 12, height = 3)

g <- ggplot(subset(tdat, mgmtdrydays == "0.1 - 0.2"), aes(x = ysincemgmt, next_x = nextysince, node = classification, next_node = nextclass, fill = factor(classification))) +
  geom_sankey(flow.alpha = 0.9, width = 0.5) + 
  facet_wrap(~mgmt, nrow = 1) + theme_linedraw() + scale_fill_manual(breaks = classes5long, values = classpal5, name = "Forest classification") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(colour = 'black')) + xlab ("Years since management") + ylab("Site-level forest classification") + 
  ggtitle("Climate impacts on forest transitions after management") + 
  geom_text(aes(label =  paste0( "n = ", size), x = Inf, y = -Inf, hjust=1, vjust=0))

ggsave("/home/ss3988/UVAFME/mgmt_results/figs/3/sitemgmt_drydays_medium_sankey.png", g, width = 12, height = 3)

g <- ggplot(subset(tdat, mgmtdrydays == "> 0.2"), aes(x = ysincemgmt, next_x = nextysince, node = classification, next_node = nextclass, fill = factor(classification))) +
  geom_sankey(flow.alpha = 0.9, width = 0.5) + 
  facet_wrap(~mgmt, nrow = 1) + theme_linedraw() + scale_fill_manual(breaks = classes5long, values = classpal5, name = "Forest classification") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(colour = 'black')) + xlab ("Years since management") + ylab("Site-level forest classification") + 
  ggtitle("Climate impacts on forest transitions after management") + 
  geom_text(aes(label =  paste0( "n = ", size), x = Inf, y = -Inf, hjust=1, vjust=0))

ggsave("/home/ss3988/UVAFME/mgmt_results/figs/3/sitemgmt_drydays_most_sankey.png", g, width = 12, height = 3)


### transparent backgrounds

g <- ggplot(subset(tdat, mgmtdrydays == "< 0.1"), aes(x = ysincemgmt, next_x = nextysince, node = classification, next_node = nextclass, fill = factor(classification))) +
  geom_sankey(flow.alpha = 0.9, width = 0.5) + 
  facet_wrap(~mgmt, nrow = 1) + scale_fill_manual(breaks = classes5long, values = classpal5, name = "Forest classification") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        strip.text = element_text(colour = 'black')) + xlab ("Years since management") + ylab("Site-level forest classification") + 
  ggtitle("Climate impacts on forest transitions after management") + 
  geom_text(aes(label =  paste0( "n = ", size), x = Inf, y = -Inf, hjust=1, vjust=0)) + transparent

ggsave("/home/ss3988/UVAFME/mgmt_results/figs/3/sitemgmt_drydays_least_sankey_trans.png", g, width = 12, height = 3)

g <- ggplot(subset(tdat, mgmtdrydays == "0.1 - 0.2"), aes(x = ysincemgmt, next_x = nextysince, node = classification, next_node = nextclass, fill = factor(classification))) +
  geom_sankey(flow.alpha = 0.9, width = 0.5) + 
  facet_wrap(~mgmt, nrow = 1) + scale_fill_manual(breaks = classes5long, values = classpal5, name = "Forest classification") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        strip.text = element_text(colour = 'black')) + xlab ("Years since management") + ylab("Site-level forest classification") + 
  ggtitle("Climate impacts on forest transitions after management") + 
  geom_text(aes(label =  paste0( "n = ", size), x = Inf, y = -Inf, hjust=1, vjust=0)) + transparent

ggsave("/home/ss3988/UVAFME/mgmt_results/figs/3/sitemgmt_drydays_medium_sankey_trans.png", g, width = 12, height = 3)

g <- ggplot(subset(tdat, mgmtdrydays == "> 0.2"), aes(x = ysincemgmt, next_x = nextysince, node = classification, next_node = nextclass, fill = factor(classification))) +
  geom_sankey(flow.alpha = 0.9, width = 0.5) + 
  facet_wrap(~mgmt, nrow = 1) +  scale_fill_manual(breaks = classes5long, values = classpal5, name = "Forest classification") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        strip.text = element_text(colour = 'black')) + xlab ("Years since management") + ylab("Site-level forest classification") + 
  ggtitle("Climate impacts on forest transitions after management") + 
  geom_text(aes(label =  paste0( "n = ", size), x = Inf, y = -Inf, hjust=1, vjust=0)) + transparent

ggsave("/home/ss3988/UVAFME/mgmt_results/figs/3/sitemgmt_drydays_most_sankey_trans.png", g, width = 12, height = 3)


#### 3.2 Forest type ~ fire probability

f = "Active_Management_Data.csv"
mgmt = data.frame(matrix(ncol=length(cols[[f]]) + 1,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen"))))
for(c in clims) {
  path = paste0(outDir, landscapeDir, "cf/", c, "/", f)
  temp = fread(path)
  temp$climscen = c
  mgmt = rbind(mgmt, temp)
}

f = "Fire_Summary.csv"
fire = data.frame(matrix(ncol=length(cols[[f]]) + 1,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen"))))
for(c in clims) {
  path = paste0(outDir, landscapeDir, "cf/", c, "/", f)
  temp = fread(path)
  temp$climscen = c
  fire = rbind(fire, temp)
}

fire = fire[!fire$year %in% c(1983, 2013),] #1983 and 2013 is double counted; initialization
fire$fprob = fire$NFires/100
dat = merge(fire, mgmt, by = c("siteID","year","climscen"))
dat$classification2 = ""
dat[dat$classification == "BSP",]$classification2 = "Black spruce"
dat[dat$classification %in% c("WSPP","WSPS"),]$classification2 = "White spruce"
dat[dat$classification %in% c("ASP","ABI","BIR"),]$classification2 = "Deciduous (aspen/birch)"
dat[dat$classification %in% c("MIX","WSB"),]$classification2 = "Mixed forest"
dat[dat$classification %in% c("REP","OTH"),]$classification2 = "Reproductive/ other"


g <- ggplot(dat, aes(x = classification2, y = fprob, fill = classification2)) + geom_boxplot() +
  theme_linedraw() + scale_fill_manual(breaks = classes5long, values = classpal5, name = "Forest classification") +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(colour = 'black')) + xlab ("Forest classification") + ylab("Fire probability") + 
  ggtitle("Annual site-level fire probability by forest type")
ggsave("/home/ss3988/UVAFME/mgmt_results/figs/3/forestclass_fireprob.png", g, width = 12, height = 8)


g <- ggplot(dat, aes(x = classification2, y = fprob, fill = classification2)) + geom_boxplot() +
  scale_fill_manual(breaks = classes5long, values = classpal5, name = "Forest classification") +
  theme(strip.text = element_text(colour = 'black')) + xlab ("Forest classification") + ylab("Fire probability") + 
  ggtitle("Annual site-level fire probability by forest type") + transparent
ggsave("/home/ss3988/UVAFME/mgmt_results/figs/3/forestclass_fireprob_trans.png", g, width = 12, height = 8)


summary(aov(NFires ~ classification2, dat))
pairwise.t.test(dat$NFires, dat$classification2)
dat %>% group_by(classification2) %>% summarize(mean(NFires))

summary(aov(I_r ~ classification2, dat))
pairwise.t.test(dat$I_r, dat$classification2)
dat %>% group_by(classification2) %>% summarize(mean(I_r, na.rm = T))


#### 3.3 Forest type ~ fire intensity


g<- ggplot(dat, aes(x = classification2, y = I_r, fill = classification2)) + geom_boxplot() +
  theme_linedraw() + scale_fill_manual(breaks = classes5long, values = classpal5, name = "Forest classification") +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(colour = 'black')) + xlab ("Forest classification") + ylab("Fire intensity (kW/m)") + 
  ggtitle("Site-level fire intensity by forest type") + 
  scale_y_continuous(trans='log10')
ggsave("/home/ss3988/UVAFME/mgmt_results/figs/3/forestclass_fireintensity.png", g, width = 12, height = 8) 

g<- ggplot(dat, aes(x = classification2, y = I_r, fill = classification2)) + geom_boxplot() +
  scale_fill_manual(breaks = classes5long, values = classpal5, name = "Forest classification") +
  theme(strip.text = element_text(colour = 'black')) + xlab ("Forest classification") + ylab("Fire intensity (kW/m)") + 
  ggtitle("Site-level fire intensity by forest type") + 
  scale_y_continuous(trans='log10') + transparent
ggsave("/home/ss3988/UVAFME/mgmt_results/figs/3/forestclass_fireintensity_trans.png", g, width = 12, height = 8) 



rm(dat)
rm(fire)
rm(mgmt)


### 4. Drivers of forest change; accelerated failure test 

# Read in data: site inputs, mgmt, clim, plot outputs
# 
# sitedat = fread(paste0(inDir, landscapeDir, "main/UVAFME2018_site.csv"))
# sitedat = sitedat[,c(1:3,6:17,37)]
# histclim = fread(paste0(inDir, landscapeDir, "main/UVAFME2018_climate.csv"))
# 
# histclim = as.data.frame(histclim %>% group_by(site) %>% mutate(prcp_win = mean(prcp_jan, prcp_feb, prcp_mar, prcp_apr, prcp_oct, prcp_nov, prcp_dec),
#                                        prcp_sum = mean(prcp_may, prcp_jun, prcp_jul, prcp_aug, prcp_sep),
#                                        tmin_win = mean(tmin_jan, tmin_feb, tmin_mar, tmin_apr, tmin_oct, tmin_nov, tmin_dec),
#                                        tmin_sum = mean(tmin_may, tmin_jun, tmin_jul, tmin_aug, tmin_sep),
#                                        tmax_win = mean(tmax_jan, tmax_feb, tmax_mar, tmax_apr, tmax_oct, tmax_nov, tmax_dec),
#                                        tmax_sum = mean(tmax_may, tmax_jun, tmax_jul, tmax_aug, tmax_sep)))
# histclim = histclim[,c(1,40:45)]
# 
# f = "Active_Management_Data.csv"
# mgmtdat = data.frame(matrix(ncol=length(cols[[f]]) + 1,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen"))))
# for(c in clims) {
#   path = paste0(outDir, landscapeDir, "cf/", c, "/", f)
#   temp = fread(path)
#   temp$climscen = c
#   mgmtdat = rbind(mgmtdat, temp)
# }
# 
# f = "Total_Plot_Values.csv"
# plotsdat = data.frame(matrix(ncol=length(cols[[f]]) + 1,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen"))))
# for(c in clims) {
#   path = paste0(outDir, landscapeDir, "cf/", c, "/", f)
#   temp = fread(path)
#   temp$climscen = c
#   plotsdat = rbind(plotsdat, temp)
# }
# 
# f = "Climate.csv"
# climdat = data.frame(matrix(ncol=length(cols[[f]]) + 1,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen"))))
# for(c in clims) {
#   path = paste0(outDir, landscapeDir, "cf/", c, "/", f)
#   temp = fread(path)
#   temp$climscen = c
#   climdat = rbind(climdat, temp)
# }
# 
# dat = merge(mgmtdat, plotsdat, by = c("siteID","year","climscen"))
# dat = merge(dat, climdat, by = c("siteID","year","climscen"))
# 
# dat = as.data.frame(dat %>% group_by(siteID, climscen) %>% mutate(prevclass = lag(classification, n = 1)))
# dat$prevclass = factor(dat$prevclass); dat$classification = factor(dat$classification)
# 
# #Need to lag annual variables because in a given year forests are classified before mortality occurs
# lagcols = c("Loreys_height","max_height","total_biomC","basal_area","total_stems", "small_stems","LAI_1","rain","pet",
#             "solar_rad","thaw_depth","organic_depth","aet", "degd","drydays")
# lagdat = dat[,c("siteID","year","climscen",lagcols)]
# lagdat = as.data.frame(lagdat %>% group_by(siteID, climscen) %>% arrange(year) %>% mutate_at(all_of(lagcols), lag, n = 1))
# colnames(lagdat)[4:18]  = paste0("prev_", colnames(lagdat)[4:18]) #this is confusing but just lazy coding...
# 
# dat = left_join(dat, lagdat, by = c("siteID","year","climscen"))
# dat = dat[,-c(6:86)]

#### 5.1 White spruce drought-driven decline; a "common" event
# before = "WSPS"
# after = c("WSB","WSPP","MIX","ABI","ASP","BIR") # white spruce declining 
# 
# wspdat = dat[dat$prevclass %in% before & dat$classification %in% after,]
# 
# #Removing changes that are transient (ie oscillating WSPS/ WSPP). Also because model below isn't designed for more than one change per site 
# wspdat = as.data.frame(wspdat %>% group_by(siteID, climscen) %>% arrange(prevclass, year) %>% slice(n()))
# wspdat$status = 1
# 
# #Now add right-censored data to wspdat
# rcensor = dat[dat$year == 2100 & dat$classification == "WSPS",]
# rcensor$status = 0
# rcensor$prevclass = "WSPS"
# wspdat = rbind(wspdat, rcensor)
# 
# #Now add environmental vars: site conditions, historic climate
# wspdat = left_join(wspdat, sitedat, by = c("siteID" = "site"))
# wspdat = left_join(wspdat, histclim[,-c(2:3)], by = c("siteID" = "site"))
# 
# wspdat$timetof= wspdat$year - wspdat$init_year
# wspdat$aspect = cos(wspdat$aspect)
# # Make training ant testing data
# set.seed(123)
# rownames(wspdat) = 1:nrow(wspdat)
# ix = sample(nrow(wspdat), round(0.7*nrow(wspdat)))
# traindat = wspdat[ix,]
# testdat = wspdat[-ix,]
# 
# y <- cbind(time=traindat$timetof, status=traindat$status)
# x <- as.matrix(traindat %>% select(-status, -timetof, -siteID, -year, -climscen, -classification, -mgmtcode, -prevclass))
# 
# cvfit <- cv.glmnet(x, y, family="cox") #first do 10-fold cross-validation to select lambda
# m <- glmnet(x, y, family="cox", lambda=cvfit$lambda.min) #plugin the optimal lambda
# cat('Cox-LASSO selected: ', paste0('X', which(m$beta!=0)), '\n')
# 
# cols = colnames(x)[c(1:4,6,8,10,12:16,18:20,22,23,25,27:32)]
# 
# #Look at which variables were selected and dropped by lasso
# cols
# colnames(x)[!colnames(x) %in% cols]
# 
# corrplot(cor(wspdat[cols]))
# wspmodel = survreg(data = traindat[c("timetof","status",cols)], formula =  Surv(timetof, status == 1, type = "right")~.) 
# summary(wspmodel)
# #Removing highly correlated and/or insignificant and/or difficult-to-observe covariates
# cols = cols[!cols %in% c("slope","aspect","prev_max_height","prev_total_biomC","prev_solar_rad","a_bd","o_bd",
#                          "prev_aet","init_year","elevation","tmin_win","latitude")]
# 
# pred <- predict(wspmodel, newdata = testdat)
# obs = testdat$timetof
# rss <- sum((pred - obs) ^ 2)  ## residual sum of squares
# tss <- sum((obs - mean(obs)) ^ 2)  ## total sum of squares
# rsq <- 1 - rss/tss
# rsq
# plot(obs, pred); abline(0,1)
# 
# modelout = as.data.frame(cbind(obs, pred))
# ggplot(modelout, aes(x = obs, y = pred, alpha = 0.6)) + geom_point() + theme_bw() + geom_abline(intercept = 0, slope = 1)
# 






sitedat = fread(paste0(inDir, landscapeDir, "main/UVAFME2018_site.csv"))
sitedat = sitedat[,c(1:3,6:17,37)]
histclim = fread(paste0(inDir, landscapeDir, "main/UVAFME2018_climate.csv"))

histclim = as.data.frame(histclim %>% group_by(site) %>% mutate(prcp_win = mean(prcp_jan, prcp_feb, prcp_mar, prcp_apr, prcp_oct, prcp_nov, prcp_dec),
                                                                prcp_sum = mean(prcp_may, prcp_jun, prcp_jul, prcp_aug, prcp_sep),
                                                                tmin_win = mean(tmin_jan, tmin_feb, tmin_mar, tmin_apr, tmin_oct, tmin_nov, tmin_dec),
                                                                tmin_sum = mean(tmin_may, tmin_jun, tmin_jul, tmin_aug, tmin_sep),
                                                                tmax_win = mean(tmax_jan, tmax_feb, tmax_mar, tmax_apr, tmax_oct, tmax_nov, tmax_dec),
                                                                tmax_sum = mean(tmax_may, tmax_jun, tmax_jul, tmax_aug, tmax_sep)))
histclim = histclim[,c(1,40:45)]

f = "Active_Management_Data.csv"
mgmtdat = data.frame(matrix(ncol=length(cols[[f]]) + 1,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen"))))
for(c in clims) {
  path = paste0(outDir, landscapeDir, "cf/", c, "/", f)
  temp = fread(path)
  temp$climscen = c
  mgmtdat = rbind(mgmtdat, temp)
}

f = "Total_Plot_Values.csv"
plotsdat = data.frame(matrix(ncol=length(cols[[f]]) + 1,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen"))))
for(c in clims) {
  path = paste0(outDir, landscapeDir, "cf/", c, "/", f)
  temp = fread(path)
  temp$climscen = c
  plotsdat = rbind(plotsdat, temp)
}

f = "Climate.csv"
climdat = data.frame(matrix(ncol=length(cols[[f]]) + 1,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen"))))
for(c in clims) {
  path = paste0(outDir, landscapeDir, "cf/", c, "/", f)
  temp = fread(path)
  temp$climscen = c
  climdat = rbind(climdat, temp)
}

dat = merge(mgmtdat, plotsdat, by = c("siteID","year","climscen"))
dat = merge(dat, climdat, by = c("siteID","year","climscen"))

dat = as.data.frame(dat %>% group_by(siteID, climscen) %>% mutate(prevclass = lag(classification, n = 1)))
dat$prevclass = factor(dat$prevclass); dat$classification = factor(dat$classification)

#Need to lag annual variables because in a given year forests are classified before mortality occurs
lagcols = c("Loreys_height","max_height","total_biomC","basal_area","total_stems", "small_stems","LAI_1","rain","pet",
            "solar_rad","thaw_depth","organic_depth","aet", "degd","drydays","flood_d")
lagdat = dat[,c("siteID","year","climscen",lagcols)]
lagdat = as.data.frame(lagdat %>% group_by(siteID, climscen) %>% arrange(year) %>% mutate_at(all_of(lagcols), lag, n = 1))
colnames(lagdat)[4:19]  = paste0("prev_", colnames(lagdat)[4:19]) #this is confusing but just lazy coding...

dat = left_join(dat, lagdat, by = c("siteID","year","climscen"))
dat = dat[,-c(6:86)]




### 5.1.1 Time-dependent variable model for wsp decline 
before = c("WSPS")
after = c("WSPP","WSB","MIX","ABI","ASP","BIR") # white spruce declining 

wspsites = unique(dat[dat$classification == before,c("siteID","climscen")])
wspdat = left_join(wspsites, dat, by = c("siteID","climscen"))
wspdat2 = dat[dat$prevclass %in% before & dat$classification %in% after,]
wspdat2 = as.data.frame(wspdat2 %>% group_by(siteID, climscen) %>% arrange(prevclass, year) %>% slice(n()))
wspdat2$status = 1

#Now add right-censored data to wspdat
rcensor = dat[dat$year == 2100 & dat$classification %in% before,]
rcensor$status = 0
rcensor$prevclass = before
wspdat2 = rbind(wspdat2, rcensor)

# Now merge back to original dataset
wspdat = left_join(wspdat, wspdat2[,c("siteID","year","climscen","status")], by = c("siteID","year","climscen"))
wspdat <- wspdat %>% mutate(status = ifelse(is.na(status),0,status))

# Censor data where wsp transitions to bsp then to decid
wspdat = wspdat[wspdat$prevclass %in% before,]
wspdat =  wspdat %>% group_by(siteID, climscen) %>% mutate(year = (1:(n())-1))

#Now add environmental vars: site conditions, historic climate
wspdat = left_join(wspdat, sitedat, by = c("siteID" = "site"))
wspdat = left_join(wspdat, histclim[,-c(2:3)], by = c("siteID" = "site"))

wspdat$aspect = cos(wspdat$aspect)
wspdat  = wspdat %>% group_by(siteID, climscen) %>% mutate(prevyear = lag(year, n = 1))
wspdat = wspdat[complete.cases(wspdat$prevyear),]
wspdat = wspdat %>% group_by(climscen, siteID) %>% mutate(id = cur_group_id()) %>% ungroup()


wspdat = as.data.frame(wspdat)
wspcols = colnames(wspdat %>% select(-siteID, -climscen, -classification, -mgmtcode, -prevclass, -init_year))
corrplot(cor(wspdat[wspcols]))
wspdat$itxt = as.factor(wspdat$itxt)
wspcols = wspcols[!wspcols %in% c("tmin_sum","o_pwp","prev_max_height","prev_total_biomC","prev_total_stems","latitude","longitude",
                                  "a_fc","o_fc","tmax_win","o_bd","o_sat","prev_pet","prev_basal_area","a_sat","itxt",
                                  "tmax_sum","prev_solar_rad","prev_aet","flood_d","prev_drydays")] #removing correlated covariates
wspmodel = coxph(formula = Surv(prevyear, year, status)~., data = wspdat[,wspcols])

wspcols = wspcols[!wspcols %in% c("slope","elevation","tmin_win","prev_rain","prev_organic_depth")] #removing very insignificant covariates
wspmodel = coxph(formula = Surv(prevyear, year, status)~., data = wspdat[,wspcols])
summary(wspmodel)


ggsurvplot(survfit(wspmodel), data = wspdat, color = "#2E9FDF", ggtheme = theme_minimal())


ggsurvplot(survfit(wspmodel), data = wspdat[wspdat$prev_Loreys_height < 15,], palette = "#2E9FDF", ggtheme = theme_minimal())
ggsave(paste0(figDir, "test.png"))
ggsurvplot(survfit(wspmodel), data = wspdat[wspdat$prev_Loreys_height >= 15,], palette = "#2E9FDF", ggtheme = theme_minimal())
ggsave(paste0(figDir, "test1.png"))

#### 5.1.2 Visualize wsp transition



#### 5.2 Black spruce transition to deciduous forest 
before = "BSP"
after = c("WSB", "WSPP","WSPS","MIX","ABI","ASP","BIR","OTH") # black spruce transitioning

bspsites = unique(dat[dat$classification == before,c("siteID","climscen")])
bspdat = left_join(bspsites, dat, by = c("siteID","climscen"))
bspdat2 = dat[dat$prevclass %in% before & dat$classification %in% after,]
bspdat2 = as.data.frame(bspdat2 %>% group_by(siteID, climscen) %>% arrange(prevclass, year) %>% slice(n()))
bspdat2$status = 1

#Now add right-censored data to bspdat
rcensor = dat[dat$year == 2100 & dat$classification == before,]
rcensor$status = 0
rcensor$prevclass = before
bspdat2 = rbind(bspdat2, rcensor)

# Now merge back to original dataset
bspdat = left_join(bspdat, bspdat2[,c("siteID","year","climscen","status")], by = c("siteID","year","climscen"))
bspdat <- bspdat %>% mutate(status = ifelse(is.na(status),0,status))
bspdat = bspdat[bspdat$prevclass == before,]
bspdat = bspdat %>% group_by(climscen, siteID) %>% mutate(id = cur_group_id()) %>% ungroup()

bspdat =  bspdat %>% group_by(siteID, climscen) %>% mutate(year = (1:(n())-1))
#Now add environmental vars: site conditions, historic climate
bspdat = left_join(bspdat, sitedat, by = c("siteID" = "site"))
bspdat = left_join(bspdat, histclim[,-c(2:3)], by = c("siteID" = "site"))

bspdat$aspect = cos(bspdat$aspect)
bspdat  = bspdat %>% group_by(siteID, climscen) %>% mutate(prevyear = lag(year, n = 1))
bspdat = bspdat[complete.cases(bspdat$prevyear),]
bspdat = as.data.frame(bspdat)
bspdat$itxt = as.factor(bspdat$itxt)
bspdat = na.omit(bspdat)

bspcols = colnames(bspdat %>% select(-siteID, -climscen, -classification, -mgmtcode, -prevclass, -init_year))
corrplot(cor((bspdat[bspcols])))
bspcols = bspcols[!bspcols %in% c("tmin_sum","o_pwp","prev_max_height","prev_total_biomC","prev_total_stems","latitude","longitude",
                                  "a_fc","o_fc","tmax_win","o_bd","o_sat","prev_pet","prev_aet","prev_basal_area","prev_Loreys_height",
                                  "prev_thaw_depth","itxt","a_sat","tmin_win","flood_d","prev_solar_rad")] #removing correlated covariates
bspmodel = coxph(formula = Surv(prevyear, year, status)~., data = bspdat[,bspcols], cluster = id)
bspcols = bspcols[!bspcols %in% c("a_bd")] #removing insignificant covariates
bspmodel = coxph(formula = Surv(prevyear, year, status)~., data = bspdat[,bspcols], cluster = id)







# # trying something new to make visualization easier
# bspcols2 = c("prev_small_stems","prev_LAI_1","prev_rain","prev_organic_depth","prev_degd","prev_drydays",
#              "elevation","slope","aspect","a_pwp","a_bd","tmax_sum")
# bspdat.cuts = surv_cutpoint(bspdat, time = "year", event = "status", variables = bspcols2)
# bspdat.cats = surv_categorize(bspdat.cuts)
# bspdat.cats$id = bspdat$id
# bspdat.cats$prevyear = bspdat$prevyear
# bspcols2 = colnames(bspdat.cats)
# bspcols2 = bspcols2[!bspcols2 %in% c("aspect","a_bd")]
# bspmodel2 = coxph(formula = Surv(prevyear, year, status)~., data = bspdat.cats[,bspcols2], cluster = id)
# bspmodel3 = coxph(formula = Surv(prevyear, year, status)~prev_small_stems + prev_LAI_1 + prev_rain + (prev_organic_depth) + prev_degd + prev_drydays +
#                     strata(elevation) + slope + a_pwp + tmax_sum, data = bspdat.cats[,bspcols2], cluster = id)
# ggsurvplot(survfit(bspmodel3), data = bspdat.cats, ggtheme = theme_minimal())
# ggforest(bspmodel2, data = bspdat.cats)
# ggsave(paste0(figDir, "4/bsp_coxph_forest.png"))
# 
# # and this for wsp:
# wspcols2 = c("prev_Loreys_height","prev_small_stems","prev_LAI_1","prev_rain","prev_thaw_depth","prev_organic_depth",
#              "prev_degd","elevation","slope","aspect","a_pwp","a_bd","tmin_win")
# wspdat.cuts = surv_cutpoint(wspdat, time = "year", event = "status", variables = wspcols2)
# wspdat.cats = surv_categorize(wspdat.cuts)
# wspdat.cats$id = wspdat$id
# wspdat.cats$prevyear = wspdat$prevyear
# wspcols2 = colnames(wspdat.cats)
# wspcols2 = wspcols2[!wspcols2 %in% c("prev_thaw_depth","prev_organic_depth")]
# newcols = c("year","status","Lorey's height (m)","Stems <5cm DBH (no./ha)",
#             "Leaf Area Index (m/m)","Precipitation (cm)","Growing degree days (C-days)",
#             "Elevation (m)","Slope (deg)","cos(aspect)",
#             "A-layer partial wilting point", "A-layer bulk density (kg/m3)",
#             "Historic winter Tmin (C)","id","prevyear")
# colnames(wspdat.cats)[colnames(wspdat.cats) %in% wspcols2] = newcols
# wspmodel2 = coxph(formula = Surv(prevyear, year, status)~., data = wspdat.cats[,newcols], cluster = id)
# ggforest(wspmodel2, data = wspdat.cats)
# ggsave(paste0(figDir, "4/wsp_coxph_forest.png"), width = 10, height = 7)
# 
# # 
# # wspmodel3 = coxph(formula = Surv(prevyear, year, status)~ prev_Loreys_height +  prev_small_stems +  strata(prev_LAI_1) +  prev_rain +
# #                   prev_degd +  elevation +  slope +  aspect +  a_pwp +  a_bd +  tmin_win, data = bspdat.cats[,wspcols2], cluster = id)
# # ggsurvplot(survfit(wspmodel3), data = wspdat.cats, ggtheme = theme_minimal())
# 
# 

# Specifying my own cutpoints because R did weird things

wspdat.cats = wspdat
wspcols = colnames(wspdat.cats %>% select(-siteID, -climscen, -classification, -mgmtcode, -prevclass, -init_year))
wspcols = wspcols[!wspcols %in% c("o_pwp","prev_max_height","prev_total_biomC","prev_total_stems","latitude","longitude",
                                  "a_fc","o_fc","tmax_win","o_bd","o_sat","prev_pet","prev_basal_area","a_sat",
                                  "prev_solar_rad","prev_aet")]  # Removing correlated rows
wspdat.cats = wspdat.cats[,wspcols]

# Select the numeric columns
num_cols <- sapply(wspdat.cats, is.numeric)
# Create histograms for the numeric columns
hist_data <- lapply(wspdat.cats[num_cols], hist)
names(hist_data)
par(mfrow = c(4,4))
for (i in c(2:10, 12:17)) {
  plot(hist_data[[i]], main = names(hist_data)[i])
}


highlow <- function(vec, threshold, a = "0", b = "1") {
  # 1 if above, 0 if below threshold
  cat <- sapply(vec, function(x) ifelse(x > threshold, a, b))
  return(cat)
}

wspdat.cats$prev_flood_d = highlow(wspdat.cats$prev_flood_d,0, "0",">0")
wspdat.cats$prev_Loreys_height = highlow(wspdat.cats$prev_Loreys_height,15, "<=15",">15")
wspdat.cats$prev_small_stems = highlow(wspdat.cats$prev_small_stems,10000, "<=10000",">10000")
wspdat.cats$prev_LAI_1 = highlow(wspdat.cats$prev_LAI_1,5.5, "<=5.5", ">5.5")
wspdat.cats$prev_rain = highlow(wspdat.cats$prev_rain,40, "<=40",">40")
wspdat.cats$prev_thaw_depth = highlow(wspdat.cats$prev_thaw_depth, 150,"<=150",">150")
wspdat.cats$prev_organic_depth = highlow(wspdat.cats$prev_organic_depth,10, "<=10",">10")
wspdat.cats$prev_degd = highlow(wspdat.cats$prev_degd,1000, "<=1000",">1000")
wspdat.cats$prev_drydays = highlow(wspdat.cats$prev_drydays,0, "0",">0")
wspdat.cats$elevation = highlow(wspdat.cats$elevation,500, "<=500",">500")
wspdat.cats$slope = highlow(wspdat.cats$slope,5, "<=5",">5")
wspdat.cats$aspect = highlow(wspdat.cats$aspect,0, "0",">0")
wspdat.cats$a_pwp = highlow(wspdat.cats$a_pwp,0.1, "<=0.1",">0.1")
wspdat.cats$a_bd = highlow(wspdat.cats$a_bd,1200, "<=1200",">1200")
wspdat.cats$tmin_win = highlow(wspdat.cats$tmin_win, -30, "<=-30",">-30")
wspdat.cats$tmax_sum = highlow(wspdat.cats$tmax_sum, 14, "<=14",">14")

wspcols2 =  c("year","prev_Loreys_height","prev_small_stems","prev_LAI_1","prev_rain","prev_thaw_depth","prev_organic_depth",
              "prev_degd","prev_drydays","prev_flood_d","status","elevation","slope","aspect","a_pwp","a_bd","itxt","tmin_win","tmax_sum",
              "prevyear","id")
newcols = c("year","Lorey's height (m)","Stems <5cm DBH (no/ha)",
            "Leaf Area Index (m/m)","Precipitation (cm)","Active layer depth (cm)","Organic layer depth (cm)",
            "Growing degree days (C-days)","Drought days (days)","Flooding days (days)","status",
            "Elevation (m)","Slope (deg)","cos(aspect)",
            "A-layer partial wilting point", "A-layer bulk density (kg/m3)","Soil texture",
            "Historic winter Tmin (C)","Historic summer Tmax (C)","prevyear","id")
colnames(wspdat.cats)[colnames(wspdat.cats) %in% wspcols2] = newcols

wspmodel3 = coxph(formula = Surv(prevyear, year, status)~., data = wspdat.cats[,newcols], cluster = id)
ggforest(wspmodel3)
ggsave(paste0(figDir, "4/wsp_coxph_forest.png"), width = 10, height = 7)



# Now for black spruce 
bspdat.cats = bspdat
bspcols = colnames(bspdat.cats %>% select(-siteID, -climscen, -classification, -mgmtcode, -prevclass, -init_year))
bspcols = bspcols[!bspcols %in% c("tmin_sum","o_pwp","prev_max_height","prev_total_biomC","prev_total_stems","latitude","longitude",
                                  "a_fc","o_fc","tmax_win","o_bd","o_sat","prev_pet","prev_aet","prev_basal_area","prev_Loreys_height",
                                  "prev_thaw_depth","a_sat","tmin_win","flood_d","prev_solar_rad")] # Removing correlated rows
bspdat.cats = bspdat.cats[,bspcols]

# Select the numeric columns
num_cols <- sapply(bspdat.cats, is.numeric)
# Create histograms for the numeric columns
hist_data <- lapply(bspdat.cats[num_cols], hist)
names(hist_data)
par(mfrow = c(4,4))
for (i in c(2:8, 11:16)) {
  plot(hist_data[[i]], main = names(hist_data)[i])
}


bspdat.cats$prev_flood_d = highlow(bspdat.cats$prev_flood_d,0, "0",">0")
bspdat.cats$prev_small_stems = highlow(bspdat.cats$prev_small_stems,20000, "<=20000",">20000")
bspdat.cats$prev_LAI_1 = highlow(bspdat.cats$prev_LAI_1,6, "<=6", ">6")
bspdat.cats$prev_rain = highlow(bspdat.cats$prev_rain,50, "<=50",">50")
bspdat.cats$prev_organic_depth = highlow(bspdat.cats$prev_organic_depth,30, "<=30",">30")
bspdat.cats$prev_degd = highlow(bspdat.cats$prev_degd,1000, "<=1000",">1000")
bspdat.cats$prev_drydays = highlow(bspdat.cats$prev_drydays,0, "0",">0")
bspdat.cats$elevation = highlow(bspdat.cats$elevation,300, "<=300",">300")
bspdat.cats$slope = highlow(bspdat.cats$slope,5, "<=5",">5")
bspdat.cats$aspect = highlow(bspdat.cats$aspect,0, "0",">0")
bspdat.cats$a_pwp = highlow(bspdat.cats$a_pwp,0.1, "<=0.1",">0.1")
bspdat.cats$a_bd = highlow(bspdat.cats$a_bd,1200, "<=1200",">1200")
bspdat.cats$tmax_sum = highlow(bspdat.cats$tmax_sum, 15, "<=15",">15")

bspcols2 =  c("year","prev_small_stems","prev_LAI_1","prev_rain","prev_organic_depth",
              "prev_degd","prev_drydays","prev_flood_d","status","id","elevation","slope","aspect","a_pwp","a_bd","itxt","tmax_sum","prevyear")

newcols = c("year","Stems <5cm DBH (no/ha)",
            "Leaf Area Index (m/m)","Precipitation (cm)","Organic layer depth (cm)",
            "Growing degree days (C-days)","Drought days (days)","Flooding days (days)","status","id",
            "Elevation (m)","Slope (deg)","cos(aspect)",
            "A-layer partial wilting point", "A-layer bulk density (kg/m3)","Soil texture",
            "Historic summer Tmax (C)","prevyear")
colnames(bspdat.cats)[colnames(bspdat.cats) %in% bspcols2] = newcols

bspmodel3 = coxph(formula = Surv(prevyear, year, status)~., data = bspdat.cats[,newcols], cluster = id)
ggforest(bspmodel3)
ggsave(paste0(figDir, "4/bsp_coxph_forest.png"), width = 10, height = 7)


#### 5.3 Export csv files with UVAFME-predicted transition year to use in QGIS
### Will map units 4C and 4D; near fairbanks, relatively intense BSP management/ WSP extraction 
# sitedat = fread(paste0(inDir, landscapeDir, "main/UVAFME2018_site.csv"))
bsptemp = bspdat #merge(bspdat, sitedat[,c("site","unit")], by.x = "siteID", by.y = "site")
wsptemp = wspdat #merge(wspdat, sitedat[,c("site","unit")], by.x = "siteID", by.y = "site")

## RCP 4.5 scenario
bsp_r45 = bsptemp[bsptemp$clim == "gcm45",]
bsp_r45 = bsp_r45[order(bsp_r45$year, decreasing = T),]
bsp_r45 = bsp_r45[order(bsp_r45$status, decreasing = T),]
bsp_r45 = bsp_r45[!duplicated(bsp_r45$siteID),]
bsp_r45 = bsp_r45[bsp_r45$unit %in% c("4C","4D"),]
bsp_r45$year = bsp_r45$year + bsp_r45$init_year + 1
fwrite(bsp_r45, paste0(figDir, "4/bsp_r45_transition.csv"))

wsp_r45 = wsptemp[wsptemp$clim == "gcm45",]
wsp_r45 = wsp_r45[order(wsp_r45$year, decreasing = T),]
wsp_r45 = wsp_r45[order(wsp_r45$status, decreasing = T),]
wsp_r45 = wsp_r45[!duplicated(wsp_r45$siteID),]
wsp_r45 = wsp_r45[wsp_r45$unit %in% c("4C","4D"),]
wsp_r45$year = wsp_r45$year + wsp_r45$init_year + 1
wsp_r45 = wsp_r45[wsp_r45$status == 1,] # none of these remain to 2100; need to filter bsp transitions
fwrite(wsp_r45, paste0(figDir, "4/wsp_r45_transition.csv"))

## RCP 8.5 scenario
bsp_r85 = bsptemp[bsptemp$clim == "gcm85",]
bsp_r85 = bsp_r85[order(bsp_r85$year, decreasing = T),]
bsp_r85 = bsp_r85[order(bsp_r85$status, decreasing = T),]
bsp_r85 = bsp_r85[!duplicated(bsp_r85$siteID),]
bsp_r85 = bsp_r85[bsp_r85$unit %in% c("4C","4D"),]
bsp_r85$year = bsp_r85$year + bsp_r85$init_year + 1
fwrite(bsp_r85, paste0(figDir, "4/bsp_r85_transition.csv"))

wsp_r85 = wsptemp[wsptemp$clim == "gcm85",]
wsp_r85 = wsp_r85[order(wsp_r85$year, decreasing = T),]
wsp_r85 = wsp_r85[order(wsp_r85$status, decreasing = T),]
wsp_r85 = wsp_r85[!duplicated(wsp_r85$siteID),]
wsp_r85 = wsp_r85[wsp_r85$unit %in% c("4C","4D"),]
wsp_r85$year = wsp_r85$year + wsp_r85$init_year + 1
wsp_r85 = wsp_r85[wsp_r85$status == 1,] # none of these remain to 2100; need to filter bsp transitions
fwrite(wsp_r85, paste0(figDir, "4/wsp_r85_transition.csv"))




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


# Supplementary figures for 500 year spinup



## New outputs are now compiled: time to generate figs and stats
### 1. Initialization vs spinup forest types (sankey plot)
dat = fread("/scratch/ss3988/UVAFME/output_data/TVSF/inits_testing/spinup/500y/Active_Management_Data.csv")
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



## New outputs are now compiled: time to generate figs and stats
### 1. Initialization vs spinup forest types (sankey plot)
dat = fread("/scratch/ss3988/UVAFME/output_data/Landsat_init/spinup/500y/Active_Management_Data.csv")
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



## Revision figures
### 1. Initialization running for 500 years at TVSF
dat = fread("/scratch/ss3988/UVAFME/output_data/TVSF/inits_testing/init/500y/Active_Management_Data.csv")
dat$classification = sapply(dat$classification, FUN = function(x) classes10hash[[toString(x)]])
years = seq(1983,2513, 5)
dat$classification = as.factor(dat$classification)
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
        axis.text = element_text(size = 15)) + biglabels + 
  ylab("Number of sites") + xlab("Year") + ggtitle("UVAFME initialized under historic climate TVSF")
ggsave(paste0(figDir, "init_500y_TVSF.png"), g, width = 10, height = 6)

### 2. Initialization running for 500 years at KLC
dat = fread("/scratch/ss3988/UVAFME/output_data/Landsat_init/spinup/init_500y/Active_Management_Data.csv")
dat$classification = sapply(dat$classification, FUN = function(x) classes11hash[[toString(x)]])
years = seq(2015, 2515, 5)
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
        axis.text = element_text(size = 15)) + biglabels + 
  ylab("Number of sites") + xlab("Year") + ggtitle("UVAFME initialized under historic climate KLC")
ggsave(paste0(figDir, "init_500y_KLC.png"), g, width = 10, height = 6)

