# #########################
# Purpose: Combine output files for Canada management results
# Author: Shelby W. Sundquist
# Date: Jan, 2023
# R version 4.1.2
# #########################
# #########################
# Input format: csv
# Output format: txt, slurm 
# #########################
# 


library(data.table)

# files to combine for later analyses
fnames = c("Total_Plot_Values.csv","Active_Management_Data.csv","Species_Data.csv", "Climate.csv",
           "Dead_Species_Data.csv","SoilDecomp.csv","Fire_Summary.csv")

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
"5 to 10","10 to 20","20 to 30","30 to 40","40 to 50","50 to 60","60 to 70",
"70 to 80","80 to 90","over 90","0 to 5 biom","5 to 10 biom","10 to 20 biom",
"20 to 30 biom","30 to 40 biom","40 to 50 biom","50 to 60 biom","60 to 70 biom",
"70 to 80 biom","80 to 90 biom","over 90 biom","degday_resp","drought_resp","shade_resp",
"perm_resp","flood_resp","nutrient_resp","max_diam","mean_diam","qmean_diam","mean_age",
"max_hgt","recr_trees","leaf_area_ind","basal_area","basal_sd","total_biomC",
"total_biomC_sd","biomC_lg","biomC_std_lg","biomC_sm","biomC_std_sm","basal_lg",
"basal_std_lg","basal_sm","basal_std_sm","dens_lg","dens_std_lg","dens_sm",
"dens_std_sm","dbh_lg","dbh_std_lg","dbh_sm","dbh_std_sm")

cols[["Climate.csv"]] = c("siteID","runID","year","rain","pet","solar_rad","thaw_depth","organic_depth",
                          "avail_n","aet","grow","pc_germ","degd","drydays","saw0_ByFC","saw0_BySAT","aow0_ByMin",
                          "wilt_days","flood_d")

cols[["Dead_Species_Data.csv"]] = c("siteID","runID","year","genus","species","0 to 5","5 to 7","7 to 10",
                                    "10 to 15","15 +","0 to 5 canker","5 to 7 canker","7 to 10 canker",
                                    "10 to 15 canker","15 + canker","0 to 5 biom","5 to 7 biom","7 to 10 biom",
                                    "10 to 15 biom","15 + biom","0 to 5 biom canker","5 to 7 biom canker",
                                    "7 to 10 biom canker","10 to 15 biom canker","15 + biom canker","degday_death",
                                    "drought_death","shade_death","perm_death","flood_death","nutrient_death",
                                    "fire_death","harv_death","wind_death","mean_diam","total_biomC","total_biomC_sd")
cols[["SoilDecomp.csv"]] = c("siteID","runID","year","odepth","odepth_sd","mdepth","moss_biom","active","OM","OM_N",
                             "alm_mining_mn","alm_mining_sd","alm_pop_mn","alm_pop_sd","aspen_mn","aspen_sd","tmin_spring",
                             "tmin_sep","flit_cornus","flit_acerfrax","flit_prunus","flit_betula","flit_queralba",
                             "flit_tsugthuj","flit_populus","flit_fagus","flit_querrubr","flit_abies","flit_picea",
                             "flit_pinus","flit_roots","flit_smboles","flit_lboles","flit_twigs","flit_smbranch",
                             "flit_lbranch","flit_WDW","flit_moss","lit_cornus","lit_acerfrax","lit_prunus",
                             "lit_betula","lit_queralba","lit_tsugthuj","lit_populus","lit_fagus","lit_querrubr",
                             "lit_abies","lit_picea","lit_pinus","lit_roots","lit_smboles","lit_lboles","lit_twigs",
                             "lit_smbranch","lit_lbranch","lit_WDW","lit_moss","avail_n")
cols[["Fire_Summary.csv"]] = c("siteID","year","max_FDI","cum_FDI","I_r","I_r_sd","NFires")
outDir = "output_data/KLC/"

cscen = c("hist","gcm45","gcm85")
groups = 20

sites = fread("input_data/KLC/main/UVAFME2018_site.csv")


for(f in fnames) {
    for(c in cscen[2:3]) {
    dat = data.frame(matrix(ncol=length(cols[[f]]),nrow=0, dimnames=list(NULL, cols[[f]])))
    for(i in 1:groups) { 
      #read all units for scenario 
      path = paste0(outDir, c, "/", i, "/", f)
      read = fread(path)
      dat = rbind(dat, read)
      print(path)
    }
    # write file with all units combined for each scenario
    path = paste0(outDir, c, "/", f)
    fwrite(dat, path)
  }
}

