#
# #########################
# Purpose: Utilities for manuscript 
# Author: Shelby W. Sundquist
# Date: Feb, 2023
# R version 4.1.2
# #########################
# #########################
# Input format: 
# Output format: 
# #########################
# 

library(data.table)
library(RColorBrewer)
library(hash)

# Directories
logDir = "logs/TVSF/"
jobDir = "jobfiles/TVSF/"
inDir = "input_data/TVSF/"
outDir = "output_data/TVSF/"
landscapeDir = "clims_testing/"
initDir = "inits_testing/"

# Vectors for looping 
amgmts = c("cf")
clims = c("hist","gcm45","gcm85")
initclims = c("","/r45","/r85")
iscen = c("init","spinup")
forests = c("bsp","wsps","decid")

# List for reading data
cols = list()
cols[["Total_Plot_Values.csv"]] = c("siteID","runID","year","gdd_death","drought_death","shade_death","perm_death","flood_death",
                                    "nutrient_death","fire_death","harvest_death","wind_death","gddresp_1","gddresp_2","gddresp_3",
                                    "gddresp_4","droughtresp_1","droughtresp_2","droughtresp_3","droughtresp_4","shaderesp_1",
                                    "shaderesp_2","shaderesp_3","shaderesp_4","permresp_1","permresp_2","permresp_3","permresp_4",
                                    "floodresp_1","floodresp_2","floodresp_3","floodresp_4","nutrientresp_1","nutrientresp_2",
                                    "nutrientresp_3","nutrientresp_4","Loreys_height","Loreys_height_sd","max_height","max_height_sd",
                                    "total_biomC","total_biomC_sd","basal_area","basal_area_sd","total_stems","total_stems_sd",
                                    "small_stems","small_stems_sd","med_stems","med_stems_sd","lg_stems","lg_stems_sd","stand_age",
                                    "stand_age_sd","LAI_1","LAI_2","LAI_3","LAI_4","LAI_5","LAI_6","LAI_7","LAI_8","LAI_9","LAI_10",
                                    "LAI_11","LAI_12")
cols[["Forestry_Data.csv"]] = c("siteID","runID","year","species","reproduction","poletimber","sawtimber","reproduction_sd",
                                "poletimber_sd","sawtimber_sd","cf_mn","cf_sd","bf_mn","bf_sd","LAI_0_2","LAI_2_4","LAI_4_6",
                                "LAI_6_8","LAI_8_10","LAI_10_12","LAI_12_14","LAI_14_16","LAI_16_18","LAI_18_20","LAI_20_22",
                                "LAI_22_24","LAI_24_26","LAI_26_28","LAI_28_30","LAI_30_32","LAI_32_34","LAI_34_36","LAI_36_38",
                                "LAI_38_40","LAI_40_60")
cols[["Active_Management_Data.csv"]] = c("siteID","year","classification","mgmtcode")
cols[["Species_Data.csv"]] = c("siteID","runID","year","genus","species","0 to 5",
                               "5 to 10","10 to 20","20 to 30","30 to 40","40 to 50",
                               "50 to 60","60 to 70","70 to 80","80 to 90","over 90",
                               "0 to 5 biom","5 to 10 biom","10 to 20 biom","20 to 30 biom",
                               "30 to 40 biom","40 to 50 biom","50 to 60 biom",
                               "60 to 70 biom","70 to 80 biom","80 to 90 biom",
                               "over 90 biom","degday_resp","drought_resp","shade_resp",
                               "perm_resp","flood_resp","nutrient_resp","max_diam",
                               "mean_diam","qmean_diam","mean_age","max_hgt","recr_trees",
                               "leaf_area_ind","basal_area","basal_sd","total_biomC",
                               "total_biomC_sd","biomC_lg","biomC_std_lg","biomC_sm",
                               "biomC_std_sm","basal_lg","basal_std_lg","basal_sm",
                               "basal_std_sm","dens_lg","dens_std_lg","dens_sm","dens_std_sm",
                               "dbh_lg","dbh_std_lg","dbh_sm","dbh_std_sm")
cols[["SoilDecomp.csv"]] = c("siteID"    ,    "runID"         ,"year"      ,    "odepth",
                             "odepth_sd"  ,   "mdepth"       , "moss_biom"   ,  "active",
                             "OM"          ,  "OM_N"          ,"alm_mining_mn", "alm_mining_sd",
                             "alm_pop_mn"   , "alm_pop_sd"   , "aspen_mn"    ,  "aspen_sd",
                             "tmin_spring"  , "tmin_sep"    ,  "flit_cornus" ,  "flit_acerfrax",
                             "flit_prunus"  , "flit_betula"  , "flit_queralba", "flit_tsugthuj",
                             "flit_populus" , "flit_fagus"  ,  "flit_querrubr", "flit_abies",
                             "flit_picea"   , "flit_pinus"   , "flit_roots"  ,  "flit_smboles",
                             "flit_lboles"  , "flit_twigs"   , "flit_smbranch" ,"flit_lbranch",
                             "flit_WDW"     , "flit_moss"    , "n_used"    ,    "avail_n")
cols[["Fire_Summary.csv"]] =  c("siteID","year","max_FDI","max_FDI_sd","I_r","I_r_sd","NFires")
cols[["Regen_Summary.csv"]] = c("siteID","year","species","regrowth","regrowth_sd")
cols[["Climate.csv"]] = c("siteID","year","runID","rain","pet","solar_rad","thaw_depth",
                          "organic_depth","avail_n","aet","grow","pc_germ","degd","summerVPD",
                          "drydays","saw0_ByFC","saw0_BySAT","aow0_ByMin","wilt_days","flood_d")
cols[["Dead_Species_Data.csv"]] = c("siteID","runID","year","genus","species","0 to 5","5 to 7","7 to 10",
                                    "10 to 15","15 +","0 to 5 canker","5 to 7 canker","7 to 10 canker",
                                    "10 to 15 canker","15 + canker","0 to 5 biom","5 to 7 biom","7 to 10 biom",
                                    "10 to 15 biom","15 + biom","0 to 5 biom canker","5 to 7 biom canker",
                                    "7 to 10 biom canker","10 to 15 biom canker","15 + biom canker","degday_death",
                                    "drought_death","shade_death","perm_death","flood_death","nutrient_death",
                                    "fire_death","harv_death","wind_death","mean_diam","total_biomC","total_biomC_sd")
# ggplot theme
specs = c("ALNUtenu","ALNUsinu","BETUneoa","PICEglau","PICEmari","POPUbals","POPUtrem","SALIscou")
specpal = c(brewer.pal(12, "Paired")[c(9,10,4,1,2,3,11,5)])

# ASP ABI BIR / MIX WSB / WSPS WSPP / BSP / REP OTH
classes5long = c("Deciduous (aspen/birch)", "Mixed forest","White spruce","Black spruce","Reproductive/ other")
classpal5 = c(brewer.pal(12, "Paired")[c(11,3,4,2,5)])

classes10 = c("ASP","ABI","BIR","BSP","MIX","WSB","WSPP","WSPS","REP","OTH")
classes10long = c("Aspen","Aspen/ birch","Birch","Black spruce","Mixed forest","White spruce/ birch","White spruce pole","White spruce saw","Reproduction","Other")
classes10hash = hash(classes10, classes10long)
classpal10 =c(brewer.pal(12, "Paired")[c(11,3,4,2,1,9,5,6,7,12)])

# climpal
gcmpal =  c(brewer.pal(11, "RdBu")[c(7,4,2)])
lccpal = c(brewer.pal(11, "RdBu")[c(7,5,4,2)])

# Visualize palettes: library(scales); show_col(palette)

TVSF_spec_theme <- function() {
  theme_bw() + 
    scale_fill_manual(breaks = specs, values = specpal) 
}

theme_small_legend <- function() {
  theme(legend.key.size = unit(0.8, 'lines'), #change legend key size
        legend.key.height = unit(.8, 'lines'), #change legend key height
        legend.key.width = unit(.8, 'lines'), #change legend key width
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=8))
}


transparent <- theme(
  strip.background = element_rect(fill = 'transparent', color = NA),
  panel.background = element_rect(fill='transparent', color = NA), #transparent panel bg
  plot.background = element_rect(fill='transparent', color = NA), #transparent plot bg
  # panel.grid.major = element_blank(), #remove major gridlines
  # panel.grid.minor = element_blank(), #remove minor gridlines
  legend.background = element_rect(fill='transparent', color = NA), #transparent legend bg
  legend.box.background = element_rect(fill='transparent') ,#transparent legend panel
  legend.key = element_rect(fill = "transparent")
  )

# Helpful functions

climnames <- function(x) {
  #x and y are vectors (e.g. one column)
  x[x == "hist"] = "Historic climate"
  x[x == "gcm45"] = "Moderate climate warming"
  x[x == "gcm85"] = "Extreme climate warming"
  return(x)
}

classes11 = c("ASP","ABI","BIR","BSP","MIX","WSB","WSPP","WSPS","REP","OTH","LPIN","APIN")
classes11long = c("Aspen","Aspen/ birch","Birch","Black spruce","Mixed forest","White spruce/ birch","White spruce pole","White spruce saw","Reproduction","Other","Lodgepole pine","Aspen/ pine")
classes11hash = hash(classes11, classes11long)
classpal11 =c(brewer.pal(12, "Paired")[c(11,3,4,2,1,9,5,6,7,12,8,10)]) #pine can be 8
classpal11 =c(brewer.pal(12, "Paired")[c(11,10,8,2,3,9,1,6,5,12,4,7)]) #pine can be 8
