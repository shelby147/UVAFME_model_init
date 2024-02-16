#Supplementary material figures

library(data.table)
library(dplyr)
library(ggplot2)
library(ggsankey)
library(reshape2)
source("scripts/analyses/KLC/Utilities.R")

spec = fread("output_data/KLC/hist/Species_Data.csv")
dead = fread("output_data/KLC/hist/Dead_Species_Data.csv")
mgmt = fread("output_data/KLC/hist/Active_Management_Data.csv")
soil = fread("output_data/KLC/hist/SoilDecomp.csv")
plots = fread("output_data/KLC/hist/Total_Plot_Values.csv")
clim = fread("output_data/KLC/hist/Climate.csv")
fire = fread("output_data/KLC/hist/Fire_Summary.csv")

spec[spec == -999] <- NA

spec45 = fread("output_data/KLC/gcm45/Species_Data.csv")
spec85 = fread("output_data/KLC/gcm85/Species_Data.csv")

spec$cscen = "hist"
spec45$cscen = "gcm45"
spec85$cscen = "gcm85"

spec = rbind(spec, spec45, spec85)
spec[spec == -999] <- NA

temp = reshape2::melt(spec[, c("siteID","year","cscen","species","degday_resp","drought_resp","shade_resp","perm_resp","flood_resp","nutrient_resp")], id.vars = c("siteID","year","cscen","species"))
temp = temp %>% group_by(cscen, year, variable, species) %>% summarize(resp = mean(value, na.rm = T), lower = quantile(value, 0.025, na.rm = T), upper = quantile(value, 0.975, na.rm = T))
temp = as.data.frame(temp)

temp$variable = factor(temp$variable, levels = c("degday_resp","drought_resp","shade_resp","perm_resp","flood_resp","nutrient_resp"), 
                       labels = c("Temperature","Drought","Shade","Permafrost","Flood","Nutrient"))

temp$cscen = factor(temp$cscen, levels = c("hist","gcm45","gcm85"), labels = c("Historical climate","SSP2-4.5", "SSP5-8.5"))


g <- ggplot(subset(temp, species == "PINUcont"), aes(x = year, y = resp, color = cscen, fill = cscen)) + facet_wrap(~variable) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) + geom_line() + labs(color = "Climate Scenario", fill = "Climate Scenario") +
  ylab("Response (0-1)") + xlab("Year") +  ggtitle("Growth-limiting factors for Kechika-Liard lodgepole pine")

ggsave(paste0(figDir, "Pine_limiting_factors.png"),g,  width= 10, height = 4)

g <- ggplot(subset(temp, species == "PICEmari"), aes(x = year, y = resp, color = cscen, fill = cscen)) + facet_wrap(~variable) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) + geom_line() + labs(color = "Climate Scenario", fill = "Climate Scenario") +
  ylab("Response (0-1)") + xlab("Year") +  ggtitle("Growth-limiting factors for Kechika-Liard black spruce")

ggsave(paste0(figDir, "Bsp_limiting_factors.png"),g,  width= 10, height = 4)


mgmt45 = fread("output_data/KLC/gcm45/Active_Management_Data.csv")
mgmt85 = fread("output_data/KLC/gcm85/Active_Management_Data.csv")

mgmt$cscen = "hist"
mgmt45$cscen = "gcm45"
mgmt85$cscen = "gcm85"

mgmt = rbind(mgmt, mgmt45, mgmt85)

temp = merge(spec, mgmt, by = c("cscen","siteID","year"))
temp = temp[temp$classification == "LPIN"]
temp = reshape2::melt(spec[, c("siteID","year","cscen","species","total_biomC")], id.vars = c("siteID","year","cscen","species"))
# temp = temp[temp$variable %in% c("degday_resp","drought_resp","shade_resp","perm_resp","flood_resp","nutrient_resp")]
temp = temp %>% group_by(cscen, year, variable, species) %>% summarize(biomC = mean(value, na.rm = T), lower = quantile(value, 0.025, na.rm = T), upper = quantile(value, 0.975, na.rm = T))
temp = as.data.frame(temp)
temp$cscen = factor(temp$cscen, levels = c("hist","gcm45","gcm85"), labels = c("Historical climate","SSP2-4.5", "SSP5-8.5"))

g <- ggplot(subset(temp[temp$species == "PINUcont",]), aes(x = year, y = biomC, color = cscen, fill = cscen)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) + geom_line() + labs(color = "Climate Scenario", fill = "Climate Scenario") +
  ylab("Biomass (tC/ha)") + xlab("Year") +  ggtitle("Kechika-Liard lodgepole pine stand biomass")

ggsave(paste0(figDir, "LPIN_biomass_timeseries.png"),g,  width= 10, height = 4)


fire45 = fread("output_data/KLC/gcm45/Fire_Summary.csv")
fire85 = fread("output_data/KLC/gcm85/Fire_Summary.csv")

fire$cscen = "hist"
fire45$cscen = "gcm45"
fire85$cscen = "gcm85"

fire = rbind(fire, fire45, fire85)
fire[fire == -999] <- NA

temp = fire %>% group_by(cscen, year) %>% summarize(fires = mean(NFires, na.rm = T))

g <- ggplot(temp, aes(x = year, y = fires, color = cscen, fill = cscen))  + geom_line()

fire$cscen = factor(fire$cscen, levels = c("hist","gcm45","gcm85"), labels = c("Historical climate","SSP2-4.5", "SSP5-8.5"))

g <- ggplot(fire, aes(year, NFires, color = cscen)) + geom_smooth() + labs(color = "Climate Scenario", fill = "Climate Scenario") +
  ylab("Site-level fires per year") + xlab("Year") +  ggtitle("Kechika-Liard annual fires per site")
ggsave(paste0(figDir, "FireFrequency.png"), g)
