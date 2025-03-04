#
# #########################
# Purpose: All figures and analyses for manuscript 2. Copying final versions for pub here so it's quick to update figs and stats if model is re-run
# Author: Shelby W. Sundquist
# Date: Jan, 2024
# R version 4.1.2
# #########################
# #########################
# Input format: txt, csv
# Output format: png, txt
# #########################
# 

# Set working directory!
path = ""
setwd(path)

library(data.table)
library(stringr)
library(sjmisc)
library(ggplot2)
library(survival)
library(ggsankey)
library(scales)
library(dplyr)
library(colorspace)
library(hash)
library(zoo)
library(glmnet)
library(survminer)
library(corrplot)
library(tidyr)
library(tidyverse)
library(gtools)
library(ggpubr)

source("scripts/analyses/Utilities.R")

figDir = "scripts/analyses/figs/"

ha_to_km2 = 0.01 # multiply
t_to_Mg = 1 #don't really need to do this, just reminding 

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

### Checking site-level outputs
for(m in cmgmts) {
  for(c in lclims) {
    fname = paste0(logDir, siteDir, "bsp/bsp_", m, "_", c, "_err.log")
    dat = readLines(fname)
    nsites = as.integer(str_extract(dat[6], "[0-9]+"))
    sitesrun = length( dat[grepl("Running for site", dat)])
    lastline = dat[length(dat)]
    if(lastline == goodlastline & nsites == sitesrun) {
      # scenario ran correctly
    } else {
      jobname = paste0("sbatch ",jobDir, "site_mgmt/bsp/bsp_", m, "_", c, "_", ".slurm")
      joberrors = c(joberrors, jobname)
      errorcodes = c(errorcodes, fname, lastline, paste0(sitesrun, " out of ", nsites, " run"))
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


## Second: Compile output files from landscape runs (site outputs will be called in their own; there are fewer and they're easier to iterate through)
# Well I finally used there/their/they're in a sentence together

# This is very time consuming; don't run it every time 

fnames = c()
# fnames = c("Total_Plot_Values.csv","Forestry_Data.csv","Active_Management_Data.csv","Climate.csv",
#            "Dead_Species_Data.csv","SoilDecomp.csv","Fire_Summary.csv","Cons_Data.csv") #"Species_Data.csv",

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
      if(f == "Cons_Data.csv") {
        # Have to fix missing column names
        colnames(dat) = c(colnames(dat)[1:17], "agw_prefire","not_burn","canopy_bd","canopy_bh","canopy_biom","R_a","rosf_active","CFB","R_final","I_final",colnames(dat)[18:19])
      }
      fwrite(dat, path)
    }
  }
}

# # combine outputs for active management; these incorporate cf scenarios for units where management doesn't occur for landscape-wide summary
print("working on active management")
for(f in fnames) {
  for(m in amgmts[2:3]) {
    for(c in clims) {
      dat = data.frame(matrix(ncol=length(cols[[f]]),nrow=0, dimnames=list(NULL, cols[[f]])))
      for(u in munits) {
        path = paste0(outDir, landscapeDir, m, "/", c, "/unit_", u, "/", f)
        read = fread(path)
        dat = rbind(dat, read)
      }
      for(u in nmunits[!nmunits %in% munits]) {
        path = paste0(outDir, landscapeDir, "cf/", c, "/unit_", u, "/", f)
        read = fread(path)
        dat = rbind(dat, read)
      }
      # write file with all units combined for each scenario
      path = paste0(outDir, landscapeDir, m, "/", c, "/", f)
            if(f == "Cons_Data.csv") {
        # Have to fix missing column names
        colnames(dat) = c(colnames(dat)[1:17], "agw_prefire","not_burn","canopy_bd","canopy_bh","canopy_biom","R_a","rosf_active","CFB","R_final","I_final",colnames(dat)[18:19])
      }
      fwrite(dat, path)
    }
  }
}

stop("Compiled files.")



# Changes to ecosystem services

## Carbon sequestration

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

mgmtdat$mgmtscen = factor(mgmtdat$mgmtscen, amgmts)
mgmtdat$climscen = factor(mgmtdat$climscen, clims)
mgmtdat$year = as.integer(mgmtdat$year)

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
plotsdat$climscen = factor(plotsdat$climscen, clims)
plotsdat$mgmtscen = factor(plotsdat$mgmtscen, amgmts)


dat = merge(plotsdat, mgmtdat, by = c("siteID","year","mgmtscen","climscen"))

dat = as.data.frame(dat)

dat$climscen = climnames2(as.character(dat$climscen))

temp = dat %>% group_by(climscen, mgmtscen, year) %>% summarize(biomass = mean(total_biomC, na.rm = T))
f1a <- ggplot(temp[temp$mgmtscen == "cf",], aes(x = year, y = biomass/ha_to_km2, color = climscen)) + geom_line(size = 1.3) + 
  xlim(2014,2100) + ylim(0,50/ha_to_km2) + theme_linedraw() + biglabels + remove_legend + ylab("Mean biomass (Mg/km2 C)") + xlab("Year") + 
  scale_color_manual(breaks = c("Historic climate","SSP2-4.5", "SSP5-8.5"), values = c("#92C5DE", "#D6604D", "#67001F"), name = "Climate scenario")
ggsave(paste0(figDir, "Fig1a.png"),f1a, width = 6, height = 4.4)

### Carbon NPV
dat = plotsdat %>% group_by(year, climscen) %>% subset(year>=2025 & mgmtscen == "cf") %>% summarize(C = mean(total_biomC)) #this is still tC/ha
dat$climscen = climnames2(as.character(dat$climscen))
## Facet b: discount rate 0%, carbon price variable

# Carbon prices from 2020-2060 from DICE 2023 model
# Money is in 2019 US dollars
Base_carbon = rep(NA, 40)
Base_carbon[1] = 6
Base_carbon[21] = 7
Base_carbon[41] = 9
Base_carbon = na.approx(Base_carbon)

Paris_carbon = rep(NA, 40)
Paris_carbon[1] = 6
Paris_carbon[21] = 47
Paris_carbon[41] = 70
Paris_carbon = na.approx(Paris_carbon)

Twodeg_carbon = rep(NA, 40)
Twodeg_carbon[1] = 6
Twodeg_carbon[6] = 18
Twodeg_carbon[21] = 155
Twodeg_carbon[41] = 306
Twodeg_carbon = na.approx(Twodeg_carbon)

# Simulating carbon value from 2025 on 
Carbonprice = dat %>% subset(year >= 2024 & year <= 2060) %>% group_by(climscen) %>% mutate(Base = C*Base_carbon[6:41])
Carbonprice = Carbonprice %>% group_by(climscen) %>% mutate(Paris = C*Paris_carbon[6:41])
Carbonprice = Carbonprice %>% group_by(climscen) %>% mutate(Twodeg = C*Twodeg_carbon[6:41])

f1b <- ggplot(Carbonprice, aes(x = year, color = climscen)) + 
  geom_line(aes(y = (Base/ha_to_km2)/1000), size = 1.3) + 
  geom_line(aes(y = (Paris/ha_to_km2)/1000), linetype = "dashed", size = 1.3) + 
  geom_line(aes(y = (Twodeg/ha_to_km2)/1000), linetype = "dotted", size = 1.3) +
  scale_y_continuous(trans='log10') + theme_linedraw() + biglabels + remove_legend +
  ylab("Carbon value (thousand 2019 USD/km2)") + xlab("Year") +#ggtitle("Carbon price scenarios and forest carbon value") + 
  scale_color_manual(breaks = c("Historic climate","SSP2-4.5", "SSP5-8.5"), values = c("#92C5DE", "#D6604D", "#67001F"), name = "Climate scenario")
ggsave(paste0(figDir, "Fig1b.png"),f1b, width = 6, height = 4)

f1b_leg = ggplot(Carbonprice, aes(x = year)) + # this figure is screwed up but the legend is okay
  geom_line(aes(y = (Base/ha_to_km2)/1000, color = climscen, linetype = "solid")) + 
  geom_line(aes(y = (Paris/ha_to_km2)/1000, color = climscen, linetype = "dashed")) + 
  geom_line(aes(y = (Twodeg/ha_to_km2)/1000, color = climscen, linetype = "dotted")) +
  scale_y_continuous(trans='log10') + theme_linedraw() + biglabels +
  ylab("Carbon value (thousand 2019 USD/km2)") + xlab("Year") +#ggtitle("Carbon price scenarios and forest carbon value") + 
  scale_linetype_manual(name = "Carbon pricing scenario",
                     labels = c("Baseline", "Paris Agreements", "2C Threshold"),
                     values = c("solid","dashed","dotted")) + 
  scale_color_manual(breaks = c("Historic climate","SSP2-4.5", "SSP5-8.5"), values = c("#92C5DE", "#D6604D", "#67001F"), name = "Climate scenario")
ggsave(paste0(figDir, "Fig1b_legend.png"),f1b_leg, width = 12, height = 8)

## Facet c: carbon price ~ paris extended, discount rate variable
Carbonprice = Carbonprice %>% group_by(climscen) %>% mutate(Discount1 = C*Paris_carbon[6:41]*(1-0.01)^(year-2025))
Carbonprice = Carbonprice %>% group_by(climscen) %>% mutate(Discount3 = C*Paris_carbon[6:41]*(1-0.03)^(year-2025))
Carbonprice = Carbonprice %>% group_by(climscen) %>% mutate(Discount5 = C*Paris_carbon[6:41]*(1-0.05)^(year-2025))


f1c <- ggplot(Carbonprice, aes(x = year, color = climscen)) + 
  geom_line(aes(y = (Discount1/ha_to_km2)/1000), size = 1.3) + 
  geom_line(aes(y = (Discount3/ha_to_km2)/1000), linetype = "dashed", size = 1.3) + 
  geom_line(aes(y = (Discount5/ha_to_km2)/1000), linetype = "dotted", size = 1.3) +
  theme_linedraw() + biglabels + remove_legend + ylab("Carbon value (2019 USD/km2") + xlab("Year") + #ggtitle("Discount rate impacts on forest carbon value") + 
  scale_color_manual(breaks = c("Historic climate","SSP2-4.5", "SSP5-8.5"), values = c("#92C5DE", "#D6604D", "#67001F"), name = "Climate scenario")
ggsave(paste0(figDir, "Fig1c.png"),f1c, width = 6, height = 4)
  
  
f1c_leg <- ggplot(Carbonprice, aes(x = year, color = climscen)) + #again, this figure is screwed up and I'm just using it for the legend
  geom_line(aes(y = Discount1/ha_to_km2, linetype = "solid")) + 
  geom_line(aes(y = Discount3/ha_to_km2, linetype = "dashed")) + 
  geom_line(aes(y = Discount5/ha_to_km2, linetype = "dotted")) +
  theme_linedraw() + biglabels + ylab("Carbon value (2019 USD/km2") + xlab("Year") + 
  scale_linetype_manual(name = "Discount rate",
                     labels = c("1%", "3%", "5%"),
                     values = c("solid","dashed","dotted")) + 
  scale_color_manual(breaks = c("Historic Climate","SSP2-4.5", "SSP5-8.5"), values = c("#92C5DE", "#D6604D", "#67001F"), name = "Climate scenario")
ggsave(paste0(figDir, "Fig1c_legend.png"),f1c_leg, width = 6, height = 4)


## Timber (WSPS) extent and yield
f = "Forestry_Data.csv"
forestdat = data.frame(matrix(ncol=length(cols[[f]]) + 1,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen"))))
m = "cf"
  for(c in clims) {
    path = paste0(outDir, landscapeDir, m, "/", c, "/", f)
    temp = fread(path)
    temp$climscen = c
    forestdat = rbind(forestdat, temp)
  }

f = "Active_Management_Data.csv"
mgmtdat = data.frame(matrix(ncol=length(cols[[f]]) + 1,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen"))))

for(c in clims) {
  path = paste0(outDir, landscapeDir, m, "/", c, "/", f)
  temp = fread(path)
  temp$climscen = c
  mgmtdat = rbind(mgmtdat, temp)
}

forestdat$climscen = factor(forestdat$climscen, clims)
dat = merge(forestdat, mgmtdat, by = c("siteID","year","climscen"))
dat$climscen = climnames2(as.character(dat$climscen))

temp = dat %>% subset(classification %in% c("WSPP","WSPS") & species == "PICEglau") %>% count(classification, climscen, year) %>% spread(classification, n, fill = 0)
temp = reshape2::melt(temp, id.vars = c("year","climscen"))
length(unique(dat$siteID)) #15819
f2a <- ggplot(temp[temp$year >= 2020,], aes(year, value/15819*100, color = climscen, linetype = variable)) + geom_line(size = 1.3) + ylim(0,15) + 
  scale_color_manual(breaks = c("Historic climate","SSP2-4.5", "SSP5-8.5"), values = c("#92C5DE", "#D6604D", "#67001F"), name = "Climate scenario") + 
  scale_linetype_manual(name = "White spruce size",
                        labels = c("Poletimber","Sawtimber"),
                        values = c(1,2)) + 
ylab("Stand type occurence (%)") + xlab("Year") + theme_linedraw() + biglabels
ggsave(paste0(figDir, "Fig2a.png"), f2a, width = 8, height = 4.5)

temp = dat %>% subset(classification %in% c("WSPP","WSPS") & species == "PICEglau") %>% group_by(year, climscen) %>% summarize(Boardfeet = mean(bf_mn))
temp$diff = NA
temp[temp$climscen == "SSP2-4.5",]$diff = temp[temp$climscen == "SSP2-4.5",]$Boardfeet / temp[temp$climscen == "Historic climate",]$Boardfeet
temp[temp$climscen == "SSP5-8.5",]$diff = temp[temp$climscen == "SSP5-8.5",]$Boardfeet / temp[temp$climscen == "Historic climate",]$Boardfeet
f2b <- ggplot(temp[temp$year >= 2020,], aes(year, diff, color = climscen)) + geom_line(size = 1.3) + ylim(0,1) + 
  scale_color_manual(breaks = c("Historic climate","SSP2-4.5", "SSP5-8.5"), values = c("#92C5DE", "#D6604D", "#67001F"), name = "Climate scenario") + 
  ylab("Volume relative to historic climate") + xlab("Year") + theme_linedraw() + biglabels + geom_hline(yintercept = 1)
ggsave(paste0(figDir, "Fig2b.png"), f2b, width = 8, height = 4.5)

temp2 = subset(temp, year >= 2025 & year <= 2045) %>% group_by(year, climscen) %>% summarize(total = sum(value))
min(temp2$total/15819)
max(temp2$total/15819)

## Habitat changes

# Coarse woody debris
m = "cf"
f = "SoilDecomp.csv"
soildat = data.frame(matrix(ncol=length(cols[[f]]) + 1,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen"))))

for(c in clims) {
  path = paste0(outDir, landscapeDir, m, "/", c, "/", f)
  temp = fread(path)
  temp$climscen = c
  soildat = rbind(soildat, temp)
}

temp = soildat
temp$br = temp$flit_smbranch + temp$flit_lbranch
temp$bl = temp$flit_smboles + temp$flit_lboles
temp = melt(temp[,c("siteID", "year", "br","bl","climscen")], id.vars = c("siteID","year","climscen"))
temp = temp %>% group_by(variable, climscen, year) %>% summarize(Vol = mean(value))
temp = dcast(temp, climscen + year ~ variable)
temp = temp[order(temp$year, temp$climscen),]
temp$dbranches = NA
temp[temp$climscen == "gcm45",]$dbranches = temp[temp$climscen == "gcm45",]$br / temp[temp$climscen == "hist",]$br
temp[temp$climscen == "gcm85",]$dbranches = temp[temp$climscen == "gcm85",]$br / temp[temp$climscen == "hist",]$br
temp$dboles = NA
temp[temp$climscen == "gcm45",]$dboles = temp[temp$climscen == "gcm45",]$bl / temp[temp$climscen == "hist",]$bl
temp[temp$climscen == "gcm85",]$dboles = temp[temp$climscen == "gcm85",]$bl / temp[temp$climscen == "hist",]$bl
temp$Branches = temp$dbranches
temp$Boles = temp$dboles
temp = melt(temp, id.vars = c("year","climscen"))
temp$climscen = climnames2(as.character(temp$climscen))

f3a <- ggplot(subset(temp, year >= 2020 & variable %in% c("Branches","Boles")), aes(year, value, linetype = variable, color = climscen)) + geom_line(size = 1.3) + 
  geom_hline(yintercept = 1) + ylab("Relative CWD weight") + xlab("Year") + theme_linedraw() + biglabels + ylim(0,1.2) +
  scale_color_manual(breaks = c("Historic climate","SSP2-4.5", "SSP5-8.5"), values = c("#92C5DE", "#D6604D", "#67001F"), name = "Climate scenario") + 
  scale_linetype_manual(name = "CWD type", labels = c("Downed branches","Dead boles"), values = c(1,2))
ggsave(paste0(figDir, "Fig3a.png"), f3a, width = 9, height = 6)

# Early seral forest

f = "Active_Management_Data.csv"
mgmtdat = data.frame(matrix(ncol=length(cols[[f]]) + 1,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen"))))

for(c in clims) {
  path = paste0(outDir, landscapeDir, m, "/", c, "/", f)
  temp = fread(path)
  temp$climscen = c
  mgmtdat = rbind(mgmtdat, temp)
}

mgmtdat$classification2 = ""
mgmtdat[mgmtdat$classification %in% c("BSP"),]$classification2 = "Black spruce"
mgmtdat[mgmtdat$classification %in% c("WSPP","WSPS","MIX","WSB"),]$classification2 = "White spruce / mixed"
mgmtdat[mgmtdat$classification %in% c("ASP","ABI","BIR"),]$classification2 = "Deciduous"
mgmtdat[mgmtdat$classification %in% c("REP"),]$classification2 = "Reproduction"
mgmtdat[mgmtdat$classification %in% c("OTH"),]$classification2 = "Other"
mgmtdat$classification2 = factor(mgmtdat$classification2, c("Other","Reproduction","Deciduous","White spruce / mixed","Black spruce"))

temp = mgmtdat %>% group_by(classification2, climscen, year) %>% summarize(count = n())
temp$climscen = climnames2(as.character(temp$climscen))
g <- ggplot(subset(temp, year >=2020 & classification2 != "Other"), aes(year, count/15819, color = climscen, linetype = classification2)) + 
  geom_line() + theme_linedraw() + biglabels + 
  ylab("Stand type occurence (%)") + xlab("Year") +
  scale_color_manual(breaks = c("Historic climate","SSP2-4.5", "SSP5-8.5"), values = c("#92C5DE", "#D6604D", "#67001F"), name = "Climate scenario") 
  
ggsave(paste0(figDir, "ForestTypeChangeLines.png"), g, width = 9, height = 6)

temp = mgmtdat %>% group_by(classification2, climscen, year) %>% summarize(count = n())
temp = temp[order(temp$climscen),]
temp1 = temp %>% subset(climscen %in% c("hist","gcm45")) %>% group_by(classification2, year) %>% mutate(rdiff1 = -1*diff(count))
temp1 = temp1 %>% subset(climscen %in% "hist") %>% group_by(classification2, year) %>% mutate(rdiff = rdiff1/count)
temp1$climscen = "gcm45"
temp2 = temp %>% subset(climscen %in% c("hist","gcm85")) %>% group_by(classification2, year) %>% mutate(rdiff1 = -1*diff(count))
temp2 = temp2 %>% subset(climscen %in% "hist") %>% group_by(classification2, year) %>% mutate(rdiff = rdiff1/count)
temp2$climscen = "gcm85"
temp3 = rbind(temp1, temp2)
temp3$climscen = climnames2(temp3$climscen)

classes5long = c("Deciduous","White spruce / mixed","Black spruce","Reproduction","Other")
classes5long = c("Other","Reproduction","Deciduous","White spruce / mixed","Black spruce")

classpal5 = c("#000000", brewer.pal(12, "Paired")[c(5,11,3,2)], "#000000") # Yellow line does not stand out on white background
classpal5 = c("#000000", brewer.pal(12, "Paired")[c(5,10,3,2)])

f3b1 <- ggplot(subset(temp, year == 2020 & climscen == "hist"), aes(year, count/15819 , fill = classification2)) + geom_col() +
  scale_fill_manual(breaks = classes5long, values = classpal5, name = "Forest type") +
  theme_linedraw() + biglabels + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  ylab("Forest type baseline proportions")
ggsave(paste0(figDir, "Fig3b1.png"), f3b1, width = 6, height = 6)

f3b2 <- ggplot(subset(temp3, classification2 != "Other" & year >= 2020), aes(year, rdiff, color = classification2, linetype = climscen)) + 
  geom_line(size = 1.3) + theme_linedraw() + biglabels + geom_hline(yintercept = 0) +
  xlab("Year") + ylab("Relative change in occurence") + 
  scale_color_manual(breaks = classes5long, values = classpal5, name = "Forest type") + 
  scale_linetype_manual(name = "Climate scenario", labels = c("SSP2-4.5", "SSP5-8.5"), values = c(1,2))
ggsave(paste0(figDir, "Fig3b2.png"), f3b2, width = 9, height = 6)

# Also need a 2020 historical baseline plot (could just be a bar with percentages in corresponding colors)

# Supplementary fig showing increase in fire mortality
f = "Total_Plot_Values.csv"
plotsdat = data.frame(matrix(ncol=length(cols[[f]]) + 1,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen"))))
m = "cf"
for(c in clims) {
  path = paste0(outDir, landscapeDir, m, "/", c, "/", f)
  temp = fread(path)
  temp$mgmtscen = m
  temp$climscen = c
  plotsdat = rbind(plotsdat, temp)
}
plotsdat$climscen = factor(plotsdat$climscen, clims)

temp = plotsdat %>% group_by(year, climscen) %>% summarize(Fmort = mean(fire_death))

f = "Fire_Conds.csv"
firedat = data.frame(matrix(ncol=length(cols[[f]]) + 1,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen"))))
m = "cf"
for(c in clims) {
  path = paste0(outDir, landscapeDir, m, "/", c, "/", f)
  temp = fread(path)
  temp$mgmtscen = m
  temp$climscen = c
  firedat = rbind(firedat, temp)
}
firedat$climscen = factor(firedat$climscen, clims)

firedat$climscen = climnames2(as.character(firedat$climscen))
temp = firedat %>% group_by(year, climscen) %>% summarize(count = n())

fsfire = ggplot(temp, aes(year, (count/15819)*100, color = climscen)) + geom_line(size = 1.3) +
  theme_linedraw() + biglabels + xlim(2025, 2100) + xlab("Year") + ylab("Annual percent of sites with fire") + ylim(0,5) +
  scale_color_manual(breaks = c("Historic climate","SSP2-4.5", "SSP5-8.5"), values = c("#92C5DE", "#D6604D", "#67001F"), name = "Climate scenario")
ggsave(paste0(figDir, "FigSFire.png"), fsfire, width = 9, height = 6)


# Supplementary fig showing what reproductive forest typically looks like 
m = "cf"
f = "Species_Data.csv"
path = paste0(outDir, landscapeDir, m, "/gcm85/", f)
specdat = fread(path)
specdat$climscen = "gcm85"

f = "Active_Management_Data.csv"
path = paste0(outDir, landscapeDir, m, "/gcm85/", f)
mgmtdat = fread(path)

dat = merge(mgmtdat, specdat, by = c("year","siteID"))
rm(specdat) # save memory
dat[dat == -999] <- NA


## Identify any regeneration failure (how long do sites stay REP?)
temp = dat[order(dat$siteID, dat$year)]
temp = temp %>% subset(species == "ALNUsinu") %>% group_by(siteID, classification) %>% mutate(ydiff = year - lag(year)) 
# if ydiff > 1, site has returned to class after being something else
# if ydiff == NA, site has never been that class before 

temp2 = temp
temp2 = temp2 %>% group_by(siteID) %>% mutate(start = ifelse(is.na(ydiff) | ydiff > 1, year, NA))
temp2 = temp2 %>% group_by(siteID) %>% mutate(end = ifelse(year == 2100, 2100, lead(start)))
# as.data.frame(temp2[temp2$siteID == unique(temp2$siteID)[1],c("year","classification","ydiff","start","end")])

starts = na.omit(temp2[,c("siteID","classification","start")])
ends = na.omit(temp2[,c("siteID","classification","end")])
temp3 = cbind(starts, ends[,"end"])
temp3$length = temp3$end - temp3$start
temp3 = temp3 %>% group_by(siteID) %>% mutate(length2 = as.factor(ifelse(end == 2100,"+","")))
# There are more "+"s than sites because some sites transitioned in 2100 and I have a headache and this isn't worth it
summary(temp3)
# hist(temp3$length)
quantile(temp3[temp3$classification == "REP",]$length, c(0.5, 0.8, 0.9, 0.95, 0.98))
# quantiles: 50% 4, 80% 10, 90% 22, 95% 38, 98% 60
# I'd say at least 80% of "reproduction stands" are true early seral and not regeneration failure in SSP 5


# Figure of what reproductive stands look like
set.seed(147)
sites = sample(unique(dat[dat$classification == "REP",]$siteID), 4)

fsrep = ggplot(subset(dat, siteID %in% sites), aes(year, total_biomC/ha_to_km2, fill = species)) + geom_area() + facet_wrap(~siteID) +
  xlim(2020, 2100) + 
  geom_vline(data = subset(temp3, siteID %in% sites & classification == "REP"), aes(xintercept = end)) +
  geom_vline(data = subset(temp3, siteID %in% sites & classification == "REP"), aes(xintercept = start)) + 
  ylab("Biomass by species (Mg/km2 C)") + xlab("Year") + 
  theme_linedraw() + biglabels + remove_legend + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(paste0(figDir, "FigSRep.png"), fsrep, width = 9, height = 6)


# Management preferences
inmgmt = fread("input_data/mgmt_testing/main/UVAFME2018_active_management.csv")
units = as.data.frame(inmgmt$unit) #  "2E"  "4C"  "4D"  "5A"  "6"   "7A"  "7B"  "10C" "14"
colnames(units) = "unit"
units$acres = c(76300, 66670, 26210, 108540, 52495, 14740, 62920, 72024, 64515)
units$ha = units$acres*0.404686

temp = melt(inmgmt, id.vars = "unit")
temp$group = c(rep("Business-as-usual", 45), rep("Adaptive silviculture: timber", 45), 
               rep("Adaptive silviculture: fuel reduction", 36))
temp$treatment = str_sub(temp$variable, start = 5L)
temp = merge(temp, units, by = "unit")
temp = temp %>% group_by(group, unit) %>% mutate(groupsum = sum(value))
temp$areapercent = round(temp$groupsum/temp$ha*100,2)

temp$treatment = factor(temp$treatment, levels = c("bspharv","bspthin","bspprune","bspshear","wsp","birch","mixed","thin","salv"),
                           labels = c("Harvest","Thin *","Prune","Shearblade","White spruce","Birch","Mixed stand","Thin","Salvage"))
# mgmtpal = brewer.pal(12,"Paired")[5:8,10,1:4]
# mgmtpal = brewer.pal(12,"Paired")[9:1]
mgmtpal = brewer.pal(12,"Paired")[1:9]

f4 = ggplot(temp, aes(x = unit, y = value*ha_to_km2, fill = treatment)) + geom_col() + ylim(0,70*ha_to_km2) + 
  facet_wrap(~group, nrow = 2) + # geom_text(aes(y = groupsum, label = paste0(areapercent,"%")), vjust = -0.5) + 
  theme_linedraw() + biglabels + ylab("Annual area treated (km2)") + xlab("Tanana Valley State Forest subunit") + 
  scale_fill_manual(values = mgmtpal) + remove_legend
ggsave(paste0(figDir, "Fig4.png"), f4, width = 9, height = 6) # manually move BAU to top left, fuels to bottom right, legend in bottom left

f4leg = ggplot(temp, aes(x = unit, y = value, fill = treatment)) + geom_col() + ylim(0,70) + 
  facet_wrap(~group, nrow = 2) + geom_text(aes(y = groupsum, label = paste0(areapercent,"%")), vjust = -0.5) + 
  theme_linedraw() + biglabels + ylab("Annual area treated (ha)") + xlab("Tanana Valley State Forest subunit") + 
  scale_fill_manual(values = mgmtpal, name = "Treatment") 
ggsave(paste0(figDir, "Fig4_legend.png"), f4leg, width = 9, height = 6)


# Management impacts at the landscape level

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

mgmtdat$mgmtscen = factor(mgmtdat$mgmtscen, amgmts)
mgmtdat$climscen = factor(mgmtdat$climscen, clims)
mgmtdat$year = as.integer(mgmtdat$year)

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
plotsdat$climscen = factor(plotsdat$climscen, clims)
plotsdat$mgmtscen = factor(plotsdat$mgmtscen, amgmts)


dat = merge(plotsdat, mgmtdat, by = c("siteID","year","mgmtscen","climscen"))
dat$climscen2 = climnames2(as.character(dat$climscen))
dat$mgmtscen2 = as.character(dat$mgmtscen)
dat$mgmtscen2 = factor(dat$mgmtscen2, levels = amgmts, labels =  c("Counterfactual","Business-as-usual","Adaptive"))


temp = dat %>% group_by(year, mgmtscen2, climscen2) %>% summarize(biom = mean(total_biomC))
fsmgmta <- ggplot(subset(temp, year > 2013), aes(year, biom/ha_to_km2, color = mgmtscen2)) + geom_line() + facet_wrap(~climscen2, nrow = 3) + ylim(0,50/ha_to_km2) + xlim(2020,2100) +
  theme_linedraw() + ylab("Mean biomass (Mg/km2 C)") + xlab("Year") + guides(color=guide_legend(title="Management Scenario")) + 
  biglabels + scale_color_brewer(palette = "Set1")
ggsave(paste0(figDir, "fsmgmta.png"), fsmgmta, height = 6, width = 8)

dat$classification2 = ""
dat[dat$classification %in% c("BSP"),]$classification2 = "Black spruce"
dat[dat$classification %in% c("WSPP","WSPS"),]$classification2 = "White spruce"
dat[dat$classification %in% c("MIX","WSB"),]$classification2 = "Mixed forest"
dat[dat$classification %in% c("ASP","ABI","BIR"),]$classification2 = "Deciduous"
dat[dat$classification %in% c("REP"),]$classification2 = "Reproduction"
dat[dat$classification %in% c("OTH"),]$classification2 = "Other"
dat$classification2 = factor(dat$classification2, c("Other","Reproduction","Deciduous","White spruce","Mixed forest","Black spruce"))

temp = dat %>% count(classification2, climscen2, mgmtscen2, year) %>% spread(classification2, n, fill = 0)
temp = reshape2::melt(temp, id.vars = c("mgmtscen2","climscen2","year"))

fsmgmtb <- ggplot(subset(temp, year > 2013 & climscen2 == "SSP2-4.5"), aes(year, value/15819*100, color = mgmtscen2)) + 
  facet_wrap(~variable, nrow = 2) + geom_line() + theme_linedraw() + xlab("Year") + ylab("Frequency on landscape (%)") + ylim(0,100) + xlim(2020,2100) + 
  guides(color="none") + biglabels + scale_color_brewer(palette = "Set1") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(paste0(figDir, "fsmgmtb.png"), fsmgmtb, height = 6, width = 9)

## How many sites are actually treated under each scenario?
# How many duplicated mgmtcodes? There are a few
temp = (mgmtdat[mgmtdat$classification == "REP" & mgmtdat$mgmtcode>0,])
temp2 = left_join(temp[,c("siteID","mgmtscen","climscen")], mgmtdat, by = c("siteID","mgmtscen","climscen"))
temp3 = temp2 %>% subset(mgmtcode > 0)

nrow(temp3 %>% subset(siteID == lag(siteID) & year == lag(year) + 1)) #90 duplicates. 75 in BAU, 15 in PRJ

mgmtdat %>% group_by(mgmtscen) %>% subset(mgmtcode > 0) %>% summarize(count = n()) # 1923 BAU, 1961 PRJ
nrow(distinct(mgmtdat[mgmtdat$mgmtscen == "bau",], siteID, year, climscen)) # 4633209
nrow(distinct(mgmtdat[mgmtdat$mgmtscen == "prj",], siteID, year, climscen)) # 4633209
# Percent area treated over simulation period
(1923-75)/4633209*100
(1961-15)/4633209*100




# Management impacts at the site level 

## Forest treatments and fuel quantity

sites = fread("input_data/site_mgmt/bsp/harv/UVAFME2018_site.csv")

f = "Fuel_Conds.csv"
fuelsdat = data.frame(matrix(ncol=length(cols[[f]]) + 2,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen","mgmtscen"))))

## This is very slow; no need to repeat it
# for(m in cmgmts) {
#   for(c in lclims) {
#     path = paste0(outDir, siteDir, "bsp/", m, "/", c, "/", f)
#     temp = fread(path)
#     # temp = temp %>% group_by(year, siteID) %>% summarize(totalFuels = max(fuel_sum), moisture = max(MEF))
#     temp = temp %>% group_by(year, siteID) %>% summarize(totalFuels = max(fuel_sum))
#     temp = temp %>% select(siteID, year, totalFuels)
#     temp$mgmtscen = m
#     temp$climscen = c
#     fuelsdat = rbind(fuelsdat, temp)
#   }
# }
# 
# fwrite(fuelsdat, paste0(outDir, siteDir, "bsp/fuelsdat.csv"))

fuelsdat = fread(paste0(outDir, siteDir, "bsp/fuelsdat.csv"))

temp = merge(fuelsdat, sites, by.x = "siteID", by.y = "site")
temp$ysincemgmt = temp$year - temp$year_management
temp2 = temp %>% group_by(ysincemgmt, mgmtscen, climscen) %>% 
  summarize(mn = mean(totalFuels), q025 = quantile(totalFuels, 0.025), q975 = quantile(totalFuels, 0.975))
temp2 = as.data.frame(temp2)
cf = subset(temp2, mgmtscen == "cf" & climscen == "0degCC")
cf = subset(cf, select = -mgmtscen)
temp2 = subset(temp2, climscen == "0degCC" & !mgmtscen %in% c("thin30","shearplant","harvplant","cf"))
temp2$mgmtscen = sapply(temp2$mgmtscen, FUN = function(x) mgmthash[[toString(x)]])

f7 <- ggplot() + 
  geom_ribbon(data = temp2, aes(x = ysincemgmt, y = mn, ymin = q025, ymax = q975, fill = mgmtscen, alpha = 0.5)) + 
  geom_line(data = temp2, aes(y = mn, x = ysincemgmt, color = mgmtscen), linewidth = 1.2) +
  geom_line(data = cf, aes(y = mn, x = ysincemgmt)) + geom_line(data = cf, aes(y = q025, x = ysincemgmt), linetype = 2) +
  geom_line(data = cf, aes(y = q975, x = ysincemgmt), linetype = 2) +
  geom_vline(xintercept = 0) +
  facet_wrap(~mgmtscen) + xlim(-1, 20) + 
  theme_linedraw() + xlab("Years since treatment") + ylab("Total fuel loading (kg/m2)") + 
  guides(fill=guide_legend(title="Treatment"), color = "none", alpha = "none") + biglabels + 
  scale_color_brewer(palette = "Set1", direction = -1) +  scale_fill_brewer(palette = "Set1", direction = -1) + remove_legend
ggsave(paste0(figDir, "Fig7.png"),f7, width = 9, height = 6)

f7leg <- ggplot() + 
  geom_ribbon(data = temp2, aes(x = ysincemgmt, y = mn, ymin = q025, ymax = q975, fill = mgmtscen, alpha = 0.5)) + 
  geom_line(data = temp2, aes(y = mn, x = ysincemgmt, color = mgmtscen), linewidth = 1.2) +
  geom_line(data = cf, aes(y = mn, x = ysincemgmt)) + geom_line(data = cf, aes(y = q025, x = ysincemgmt), linetype = 2) +
  geom_line(data = cf, aes(y = q975, x = ysincemgmt), linetype = 2) +
  geom_vline(xintercept = 0) +
  facet_wrap(~mgmtscen) + xlim(-1, 20) + 
  theme_linedraw() + xlab("Years since treatment") + ylab("Total fuel loading (kg/m^-2)") + 
  guides(fill=guide_legend(title="Treatment"), color = "none", alpha = "none") + biglabels + 
  scale_fill_brewer(palette = "Set1", direction = -1) + scale_color_brewer(palette = "Set1", direction = -1)
ggsave(paste0(figDir, "Fig7legend.png"),f7leg, width = 9, height = 6)

## Forest treatments and forest change (deciduous fraction)

#### Deciduous fraction is DF by basal area

### intermediate data product saved below, because it takes a while to read all the dependent files
# 
# sites = fread("input_data/site_mgmt/bsp/harv/UVAFME2018_site.csv")
# 
# f = "Climate.csv"
# climdat = data.frame(matrix(ncol=length(cols[[f]]) + 2,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen","mgmtscen"))))
# for(m in cmgmts[c(1,2,4,6,8)]) {
#   for(c in lclims) {
#     path = paste0(outDir, siteDir, "bsp/", m, "/", c, "/", f)
#     temp = fread(path)
#     temp$mgmtscen = m
#     temp$climscen = c
#     climdat = rbind(climdat, temp)
#   }
# }
# 
# f = "Active_Management_Data.csv"
# mgmtdat = data.frame(matrix(ncol=length(cols[[f]]) + 2,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen","mgmtscen"))))
# for(m in cmgmts[c(1,2,4,6,8)]) {
#   for(c in lclims) {
#     path = paste0(outDir, siteDir, "bsp/", m, "/", c, "/", f)
#     temp = fread(path)
#     temp$mgmtscen = m
#     temp$climscen = c
#     mgmtdat = rbind(mgmtdat, temp)
#   }
# }
# 
# dat = merge(mgmtdat, climdat, by = c("siteID","year","mgmtscen","climscen"))
# dat = merge(dat, sites[,c("site","year_management","latitude","longitude","aspect","elevation")], by.x = "siteID", by.y = "site")
# dat$ysincemgmt = dat$year - dat$year_management
# 
# f = "Species_Data.csv"
# specdat = data.frame(matrix(ncol=length(cols[[f]]) + 2,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen","mgmtscen"))))
# for(m in cmgmts[c(1,2,4,6,8)]) {
#   for(c in lclims) {
#     path = paste0(outDir, siteDir, "bsp/", m, "/", c, "/", f)
#     temp = fread(path)
#     temp$mgmtscen = m
#     temp$climscen = c
#     specdat = rbind(specdat, temp)
#   }
# }
# 
# f = "Total_Plot_Values.csv"
# plotsdat = data.frame(matrix(ncol=length(cols[[f]]) + 2,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen","mgmtscen"))))
# for(m in cmgmts[c(1,2,4,6,8)]) {
#   for(c in lclims) {
#     path = paste0(outDir, siteDir, "bsp/", m, "/", c, "/", f)
#     temp = fread(path)
#     temp$mgmtscen = m
#     temp$climscen = c
#     plotsdat = rbind(plotsdat, temp)
#   }
# }
# 
# specdat[specdat == -999] <- NA
# specdat$decid = 0
# specdat[specdat$species %in% c("ALNUsuni","ALNUtenu","POPUbals","POPUtrem","BETUneoa","SALIscou")]$decid = 1
# DF = specdat %>% group_by(year, siteID, decid, mgmtscen, climscen) %>% summarize(sum(basal_area))

# DF = merge(DF, plotsdat[,c("siteID","year","basal_area","mgmtscen","climscen")], by = c("siteID","year","mgmtscen","climscen"))
# DF = DF[DF$decid == 1,]
# DF$DF = DF$`sum(basal_area)`/DF$basal_area

# dat = merge(as.data.frame(dat), DF[,c("siteID","year","mgmtscen","climscen","DF")], by = c("siteID","year","mgmtscen","climscen"))

# specdat$decid = 0
# specdat[specdat$species %in% c("POPUbals","POPUtrem","BETUneoa")]$decid = 1
# DF = specdat %>% group_by(year, siteID, decid, mgmtscen, climscen) %>% summarize(sum(basal_area))

# DF = merge(DF, plotsdat[,c("siteID","year","basal_area","mgmtscen","climscen")], by = c("siteID","year","mgmtscen","climscen"))
# DF = DF[DF$decid == 1,]
# DF$DF2 = DF$`sum(basal_area)`/DF$basal_area

# dat = merge(as.data.frame(dat), DF[,c("siteID","year","mgmtscen","climscen","DF2")], by = c("siteID","year","mgmtscen","climscen"))
# fwrite(dat, paste0(outDir, siteDir, "bsp/DeciduousFraction.csv"))


dat = fread(paste0(outDir, siteDir, "bsp/DeciduousFraction.csv"))


# site-level deciduous fraction differences
temp2 = reshape2::dcast(dat[,c("siteID","ysincemgmt","mgmtscen","climscen","DF")],  siteID + climscen + ysincemgmt ~ mgmtscen)
temp2$harvdiff = temp2$harv - temp2$cf
temp2$sheardiff = temp2$shear - temp2$cf
temp2$prunediff = temp2$prune - temp2$cf
temp2$thindiff = temp2$thin - temp2$cf

temp3 = reshape2::melt(temp2[,-c(4:8)], id.vars = c("ysincemgmt","climscen","siteID"))
temp3 = temp3 %>% group_by(ysincemgmt, climscen, variable) %>% 
  summarize(diff = mean(value, na.rm = T), q025 = quantile(value, 0.25, na.rm = T), q975 = quantile(value, 0.75, na.rm = T))
temp3 = temp3[temp3$ysincemgmt != 0,]

fsdf <- ggplot(temp3[temp3$climscen == "0degCC",], aes(x = ysincemgmt, y = diff, color = variable))  + 
  xlim(-1, 40) + ylim(0,0.5) + geom_line(size = 1.3) + ylab("Within-site deciduous fraction change") + xlab("Years since treatment") + 
  theme_linedraw() + biglabels +   scale_color_brewer(palette = "Set1", direction = -1) + remove_legend + geom_hline(yintercept = 0)
ggsave(paste0(figDir, "FigSDF.png"), fsdf, width = 9, height = 6)


# DF change when not considering shrubs
temp2 = reshape2::dcast(dat[,c("siteID","ysincemgmt","mgmtscen","climscen","DF2")],  siteID + climscen + ysincemgmt ~ mgmtscen)
temp2$harvdiff = temp2$harv - temp2$cf
temp2$sheardiff = temp2$shear - temp2$cf
temp2$prunediff = temp2$prune - temp2$cf
temp2$thindiff = temp2$thin - temp2$cf

temp3 = reshape2::melt(temp2[,-c(4:8)], id.vars = c("ysincemgmt","climscen","siteID"))
temp3 = temp3 %>% group_by(ysincemgmt, climscen, variable) %>% 
  summarize(diff = mean(value, na.rm = T), q025 = quantile(value, 0.25, na.rm = T), q975 = quantile(value, 0.75, na.rm = T))
temp3 = temp3[temp3$ysincemgmt != 0,]

fsdfb <- ggplot(temp3[temp3$climscen == "0degCC",], aes(x = ysincemgmt, y = diff, color = variable))  + 
  # geom_ribbon(aes(ymin = q025, ymax = q975, alpha = 0.1, fill = NA, color = variable)) + 
  xlim(-1, 40) + geom_line(size = 1.3) + ylab("Within-site deciduous fraction change") + xlab("Years since treatment") + 
  theme_linedraw() + biglabels +   scale_color_brewer(palette = "Set1", direction = -1) + remove_legend + geom_hline(yintercept = 0)
ggsave(paste0(figDir, "FigSDFb.png"), fsdfb, width = 9, height = 6)


## Forest type and fire frequency/ intensity

f = "Active_Management_Data.csv"
mgmt = data.frame(matrix(ncol=length(cols[[f]]) + 1,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen"))))
for(c in clims) {
  path = paste0(outDir, landscapeDir, "cf/", c, "/", f)
  temp = fread(path)
  temp$climscen = c
  mgmt = rbind(mgmt, temp)
}

# Final fire intensity
f = "Cons_Data.csv" # the fires that actually happen
firedat = data.frame(matrix(ncol=length(cols[[f]]) + 1,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen"))))
m = "cf"
for(c in clims) {
  path = paste0(outDir, landscapeDir, m, "/", c, "/", f)
  temp = fread(path)
  temp$mgmtscen = m
  temp$climscen = c
  firedat = rbind(firedat, temp)
}

dat = merge(firedat, mgmt, by = c("siteID","year","climscen"))

temp = dat %>% group_by(climscen, classification) %>% summarize(Ifinal_mn = mean(I_final), Ifinal_975 = quantile(I_final, 0.975), Ifinal_025 = quantile(I_final, 0.025))
temp2 = mgmt %>% group_by(classification, climscen) %>% summarize(count = n())
temp2$area_freq = temp2$count/nrow(mgmt)*3 #multiply by 3 for each climate scenario. Climscens are represented equally in mgmt

temp3 = merge(temp, temp2, by = c("climscen","classification"))

# Shoulda done this before the merge, oh well
dat$classification2 = ""
dat[dat$classification %in% c("BSP"),]$classification2 = "Black spruce"
dat[dat$classification %in% c("WSPP","WSPS"),]$classification2 = "White spruce"
dat[dat$classification %in% c("MIX","WSB"),]$classification2 = "Mixed forest"
dat[dat$classification %in% c("ASP","ABI","BIR"),]$classification2 = "Deciduous"
dat[dat$classification %in% c("REP"),]$classification2 = "Reproduction"
dat[dat$classification %in% c("OTH"),]$classification2 = "Other"
dat$classification2 = factor(dat$classification2, c("Other","Reproduction","Deciduous","White spruce","Mixed forest","Black spruce"))

mgmt$classification2 = ""
mgmt[mgmt$classification %in% c("BSP"),]$classification2 = "Black spruce"
mgmt[mgmt$classification %in% c("WSPP","WSPS"),]$classification2 = "White spruce"
mgmt[mgmt$classification %in% c("MIX","WSB"),]$classification2 = "Mixed forest"
mgmt[mgmt$classification %in% c("ASP","ABI","BIR"),]$classification2 = "Deciduous"
mgmt[mgmt$classification %in% c("REP"),]$classification2 = "Reproduction"
mgmt[mgmt$classification %in% c("OTH"),]$classification2 = "Other"
mgmt$classification2 = factor(mgmt$classification2, c("Other","Reproduction","Deciduous","White spruce","Mixed forest","Black spruce"))

temp = dat %>% subset(climscen == "gcm45" & year >=25 & year <=2050) %>% group_by(classification2) %>% summarize(Ifinal_mn = mean(I_final), Ifinal_med = median(I_final),
                                                                                                                 Ifinal_975 = quantile(I_final, 0.975), Ifinal_025 = quantile(I_final, 0.025),
                                                                                                                 Ifinal_25 = quantile(I_final, 0.25), Ifinal_75 = quantile(I_final, 0.75))                                                 
temp2 = mgmt %>% subset(climscen == "gcm45" & year >=2025 & year <=2050) %>% group_by(classification2) %>% summarize(count = n())
temp2$area_freq = temp2$count/nrow(mgmt[mgmt$year >=2025 & mgmt$year <= 2050 & mgmt$climscen == "gcm45",]) #multiply by 3 for each climate scenario. Climscens are represented equally in mgmt

temp3 = merge(temp, temp2, by = c("classification2"))

classpal5 = c(brewer.pal(12, "Paired")[c(5,10,3,2,8)])
classes5long = c("Reproduction","Deciduous","Mixed forest","Black spruce","White spruce")

temp4 = merge(dat[dat$climscen == "gcm45" & dat$year >=2025 & dat$year <=2050,], temp2, by= "classification2")

f6v2 <- ggplot(subset(temp4, climscen == "gcm45" & year >=2025 & year <=2050 & classification2 != "Other"), 
               aes(x = area_freq*100, y = I_final/1000, color = classification2)) + 
  geom_point(position = "jitter", alpha = 0.25) + xlim(5,35) +
  theme_linedraw() + biglabels + ylab("Fireline intensity (MW/m)") + xlab("Percent area by forest type") + 
  scale_color_manual(breaks = classes5long, values = classpal5, name = "Forest type") + remove_legend
# geom_violin()
ggsave(paste0(figDir, "Fig6v2.png"), f6v2, height = 5, width = 7)


f6legv2 <- ggplot(subset(temp4, climscen == "gcm45" & year >=2025 & year <=2050 & classification2 != "Other"), 
               aes(x = area_freq*100, y = I_final/1000, color = classification2)) + 
  geom_point(position = "jitter", size = 3) + xlim(0,35) +
  theme_linedraw() + biglabels + ylab("Fire intensity (MW/m)") + xlab("Percent area by forest type") + 
  scale_color_manual(breaks = classes5long, values = classpal5, name = "Forest type") 
# geom_violin()
ggsave(paste0(figDir, "Fig6v2_legend.png"), f6legv2, height = 6, width = 9)


# Explore bsp fire intensity bimodal distribution
# Identify low-point in bsp bimodal distribution
d <- subset(temp4, climscen == "gcm45" & year >=2025 & year <=2050 & classification2 == "Black spruce")
d = density(d$I_final)
optimize(approxfun(d),interval=c(0,40000))$minimum #11515.59

g <- ggplot(subset(temp4, climscen == "gcm45" & year >=2025 & year <=2050 & classification2 == "Black spruce"), 
            aes(x = I_final)) + geom_histogram(binwidth = 5000) + geom_vline(xintercept = 11515.59)
ggsave(paste0(figDir, "test3.png"), g, height = 6, width = 9) # need to check this figure and make sure low point line is in the right place

tsfire = temp3 %>% summarize_if(is.numeric, mean)
tsfire$classification2 = "All"; tsfire$area_freq = sum(temp3$area_freq); tsfire$count = sum(temp3$count)
tsfire = rbind(temp3, tsfire)

temp5 = temp4 %>% subset(classification2 == "Black spruce") %>% mutate(bspgroup = ifelse(I_final > 11515.59, "High","Low"))
summary(as.factor(temp5$bspgroup)) # high: 339, Low: 2393
temp6 = temp5 %>% group_by(bspgroup) %>% summarize(Ifinal_mn = mean(I_final), Ifinal_med = median(I_final), 
                                                 Ifinal_975 = quantile(I_final, 0.975), Ifinal_025 = quantile(I_final, 0.025),
                                                 Ifinal_25 = quantile(I_final, 0.25), Ifinal_75 = quantile(I_final, 0.75))
                                                                                                                 
temp6$count = NA
temp6$area_freq = NA
temp6$classification2 = c("High intensity","Low intensity")

tsfire = rbind(tsfire, temp6[,-1])
tsfire = tsfire[c(7,1,9,8,6,3,2,5,4),]
fwrite(tsfire, paste0(figDir, "tsfire_v2.csv"))


# Need to go get these other variables again
f = "Fire_Conds.csv" # the fires that actually happen
firedat2 = data.frame(matrix(ncol=length(cols[[f]]) + 1,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen"))))
m = "cf"
for(c in clims) {
  path = paste0(outDir, landscapeDir, m, "/", c, "/", f)
  temp = fread(path)
  temp$mgmtscen = m
  temp$climscen = c
  firedat2 = rbind(firedat2, temp)
}

dat2 = merge(firedat2, mgmt, by = c("siteID","year","climscen"))
dat2$classification2 = ""
dat2[dat2$classification %in% c("BSP"),]$classification2 = "Black spruce"
dat2[dat2$classification %in% c("WSPP","WSPS"),]$classification2 = "White spruce"
dat2[dat2$classification %in% c("MIX","WSB"),]$classification2 = "Mixed forest"
dat2[dat2$classification %in% c("ASP","ABI","BIR"),]$classification2 = "Deciduous"
dat2[dat2$classification %in% c("REP"),]$classification2 = "Reproduction"
dat2[dat2$classification %in% c("OTH"),]$classification2 = "Other"
dat2$classification2 = factor(dat2$classification2, c("Other","Reproduction","Deciduous","White spruce","Mixed forest","Black spruce"))

temp4_2 = merge(dat2[dat2$climscen == "gcm45" & dat2$year >=2025 & dat2$year <=2050,], temp2, by= "classification2")
temp5_2 = merge(temp4_2, temp5[,c("siteID","year","climscen","I_final","bspgroup")], all.y = T, by = c("year","climscen","siteID"))

bspfires = melt(temp5_2[,-c("classification","mgmtscen","climscen","classification2")], id.vars = c("bspgroup"))
vars = c("DRI","day","MEF","fuel_moss","fuel_mosslive","fuel_dec","fuel_con","fuel_bole","fuel_smbr")

fsfuel <- ggplot(subset(bspfires, variable %in% vars), aes(fill = bspgroup, x = value)) + 
  facet_wrap(~variable,scales = "free") + geom_histogram(position="identity", alpha = 0.5) + 
  theme_linedraw() + biglabels + scale_fill_manual(name = "Fire Intensity", breaks = c("Low","High"), values = c("#1F78B4","#E31A1C")) + 
  remove_legend
ggsave(paste0(figDir, "fsfuel_v2.png"), fsfuel, height = 8, width = 12)


fsfuel_leg <- ggplot(subset(bspfires, variable %in% vars), aes(fill = bspgroup, x = value)) + 
  facet_wrap(~variable,scales = "free") + geom_histogram(position="identity", alpha = 0.5) + 
  theme_linedraw() + biglabels + scale_fill_manual(name = "Fireline Intensity", breaks = c("Low","High"), values = c("#1F78B4","#E31A1C"))
ggsave(paste0(figDir, "fsfuel_v2_leg.png"), fsfuel_leg, height = 6, width = 9)


# White spruce replanting success

# Supplementary fig showing what reproductive forest typically looks like
m = "bau"
f = "Species_Data.csv"
path = paste0(outDir, landscapeDir, m, "/hist/", f)
specdat = fread(path)

f = "Active_Management_Data.csv"
path = paste0(outDir, landscapeDir, m, "/hist/", f)
mgmtdat = fread(path)

dat = merge(mgmtdat, specdat, by = c("year","siteID"))
rm(specdat) # save memory
dat[dat == -999] <- NA

sites = unique(dat[dat$mgmtcode == 1 & dat$classification == "WSPS" & year <=2040]$siteID)
set.seed(147*2)
s = sample(sites, 6)

temp = dat %>% subset(siteID %in% s) %>% group_by(siteID) %>% mutate(mgmtyear = ifelse(mgmtcode == 1, year, NA))
temp = temp %>% group_by(siteID) %>% mutate(wspagain = ifelse(year>mgmtyear & classification %in% c("WSPS","WSPP","WSB"),year,NA))
temp$wspagain = as.numeric(temp$wspagain)

fsreplant = ggplot(temp, aes(year, total_biomC/ha_to_km2, fill = species)) + facet_wrap(~siteID) + geom_area() + geom_vline(aes(xintercept = mgmtyear)) + xlim(2014, 2100) + 
 theme_linedraw() + biglabels + xlab("Year") + ylab("Biomass (Mg/km2 C) by species") + remove_legend +   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(paste0(figDir, "SReplanting.png"), fsreplant, height = 6, width = 9)

# How many sites are wsps 40 years after replanting in each climate scenario?
f = "Active_Management_Data.csv"
mgmtdat = data.frame(matrix(ncol=length(cols[[f]]) + 2,nrow=0, dimnames=list(NULL, c(cols[[f]], "climscen","mgmtscen"))))

for(m in amgmts[2:3]) {
  for(c in clims) {
    path = paste0(outDir, landscapeDir, m, "/", c, "/", f)
    temp = fread(path)
    temp$mgmtscen = m
    temp$climscen = c
    mgmtdat = rbind(mgmtdat, temp)
  }
}

dat = mgmtdat[0,]
for(m in amgmts[2:3]) {
  for(c in clims) {
  sites = unique(mgmtdat[mgmtdat$mgmtcode == 1 & mgmtdat$classification == "WSPS" & mgmtdat$climscen == c & mgmtdat$mgmtscen == m,]$siteID)
  temp = mgmtdat %>% subset(climscen == c & mgmtscen == m & siteID %in% sites)
  dat = rbind(dat, temp)
  }
}

temp = dat %>% group_by(climscen, mgmtscen, siteID) %>% mutate(mgmtyear = ifelse(mgmtcode == 1, year, NA))
temp = temp %>% group_by(climscen, mgmtscen, siteID) %>% mutate(mgmtyear = ifelse(is.na(mgmtyear), max(mgmtyear, na.rm = T), mgmtyear))
# temp = temp %>% subset(mgmtyear <= 2060) # giving 40 years for stand to return to wspp/ wsps/ wsb
# temp2 = temp %>% subset(year > mgmtyear & classification %in% c("WSPS","WSPP","WSB"))
temp3 = temp %>% subset(year > mgmtyear & year <= mgmtyear + 40 & classification %in% c("WSPS","WSPP","WSB"))
nrow(distinct(temp3[,c("siteID","climscen","mgmtscen")]))
temp3_1 = temp %>% subset(temp$mgmtyear <=2100-40)
nrow(distinct(temp3_1[,c("siteID","climscen","mgmtscen")]))

nrow(distinct(temp3[temp3$climscen == "hist",c("siteID","climscen","mgmtscen")]))
nrow(distinct(temp3_1[temp3_1$climscen == "hist",c("siteID","climscen","mgmtscen")]))

nrow(distinct(temp3[temp3$climscen == "gcm85",c("siteID","climscen","mgmtscen")]))
nrow(distinct(temp3_1[temp3_1$climscen == "gcm85",c("siteID","climscen","mgmtscen")]))

# how many within 60 years?
temp4 = temp %>% subset(year > mgmtyear & year <= mgmtyear + 60 & classification %in% c("WSPS","WSPP","WSB"))
nrow(distinct(temp4[,c("siteID","climscen","mgmtscen")]))
temp4_1 = temp %>% subset(temp$mgmtyear <=2100-60)
nrow(distinct(temp4_1[,c("siteID","climscen","mgmtscen")]))


nrow(distinct(temp4[temp4$climscen == "hist",c("siteID","climscen","mgmtscen")]))
nrow(distinct(temp4_1[temp4_1$climscen == "hist",c("siteID","climscen","mgmtscen")]))

nrow(distinct(temp4[temp4$climscen == "gcm85",c("siteID","climscen","mgmtscen")]))
nrow(distinct(temp4_1[temp4_1$climscen == "gcm85",c("siteID","climscen","mgmtscen")]))
