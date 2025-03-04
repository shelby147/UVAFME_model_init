
#
# #########################
# Purpose: Set up input files to run black spruce fuels treatments test in UVAFME
# Author: Shelby W. Sundquist
# Date: Jan, 2024
# R version 4.1.2
# #########################
# #########################
# Input format: csv
# Output format: csv, slurm
# #########################
# 

# Set working directory!
setwd("")
library(data.table)
library(dplyr)
library(lubridate)


# New site files
# bsp treatments: shearblade, prune, thin30-60 percent, clearcut, counterfactual
outDir = "input_data/site_mgmt/bsp/"
slist = fread(paste0(outDir, "UVAFME2018_sitelist.csv"))
sites = as.data.frame(fread("input_data/mgmt_testing/main/UVAFME2018_site.csv"))
sites = sites[sites$site %in% slist$site,]

sites$management = 0
sites$sel_cutting = 0
sites$planting = "nothingg"
sites$viability = 0
sites$rotation_time = 0
sites$dbh_thresh = 0
sites$thinning = 0
sites$year_management = 0
sites$thin_perc = 0.0
sites$shearblading = 0
sites$pruning = 0
sites$init_year = 2013

sitescf = sites

set.seed(147)
sites$year_management = round(runif(500, 2014, 2080))
sites$management = 1

sitesharv = sites 
sitesharv$sel_cutting = 1
sitesharv$dbh_thresh = 0
sitesharv$rotation_time = 200 # just setting this high, not interested

sitesharvplant = sitesharv
sitesharvplant$planting = "PICEglau"
sitesharvplant$viability = 1

sitesthin30 = sites
sitesthin30$thinning = 1
sitesthin30$thin_perc = .30

sitesthin60 = sitesthin30
sitesthin60$thin_perc = .60

sitesshear = sites
sitesshear$shearblading = 1
sitesshearplant = sitesshear
sitesshearplant$planting = "PICEglau"
sitesshearplant$viability = 1

sitesprune = sites
sitesprune$pruning = 1

dir.create(file.path(paste0(outDir, 'cf')), showWarnings = FALSE)
dir.create(file.path(paste0(outDir,'prune')), showWarnings = FALSE)
dir.create(file.path(paste0(outDir,'shear')), showWarnings = FALSE)
dir.create(file.path(paste0(outDir,'thin30')), showWarnings = FALSE)
dir.create(file.path(paste0(outDir,'thin60')), showWarnings = FALSE)
dir.create(file.path(paste0(outDir,'harv')), showWarnings = FALSE)
dir.create(file.path(paste0(outDir,'harvplant')), showWarnings = FALSE)
dir.create(file.path(paste0(outDir,'shearplant')), showWarnings = FALSE)

fwrite(sites, paste0(outDir,'cf/UVAFME2018_site.csv'))
fwrite(sitesharv, paste0(outDir,'harv/UVAFME2018_site.csv'))
fwrite(sitesharvplant, paste0(outDir,'harvplant/UVAFME2018_site.csv'))
fwrite(sitesthin30, paste0(outDir,'thin30/UVAFME2018_site.csv'))
fwrite(sitesthin60, paste0(outDir,'thin60/UVAFME2018_site.csv'))
fwrite(sitesshear, paste0(outDir,'shear/UVAFME2018_site.csv'))
fwrite(sitesshearplant, paste0(outDir,'shearplant/UVAFME2018_site.csv'))
fwrite(sitesprune, paste0(outDir,'prune/UVAFME2018_site.csv'))


# New filelists
l1 = "&filenames"
l2 = "input_directory='input_data/mgmt_testing/main/'"
l3 = "output_directory='output_data/site_mgmt/" #/_mgmt_/_x_degCC
l4 = "climate_directory='input_data/mgmt_testing/main/'"
l5 = "site_directory='input_data/site_mgmt/bsp/" #_mgmt_
l6 = "sitelist_directory='input_data/site_mgmt/"
l7 = "GCM_directory='input_data/mgmt_testing/main/'"
l8 = "rt_directory='input_data/site_mgmt/bsp/rtdir/" #_x_degCC
l9 = "speclist_directory='input_data/mgmt_testing/main/'"
l10 = "/"

CC = c(0,1,2,4) # Exoeriment with warming/ drying forcing
CC = 0 #historical climate only 
mgmts = c("cf","harv","harvplant","thin30","thin60","shear","shearplant","prune")

f = "bsp"
for(m in mgmts) {
  for(deg in CC) {
    
    dir = paste0(f, "/", m,"/",deg,"degCC")
    l3a = paste0(l3, dir,"'")
    # Create output directory
    # dir.create(file.path(paste0("output_data/site_mgmt/", dir)), showWarnings = FALSE)
    
    l5a = paste0(l5, m,"'")
    l6a = paste0(l6, f,"'")
    l8a = paste0(l8, deg, "degCC'")
    
    
    flist = c(l1, l2, l3a, l4, l5a, l6a, l7, l8a, l9, l10)
    writeLines(flist, paste0("file_lists/site_mgmt/bsp/file_list_",f,"_",m,"_",deg,"degCC.txt"))
  }
}


# New batch scripts with increased time allowance

jobdir = "jobfiles/site_mgmt/bsp/"
model = "UVAFME_init"
joblist = c()

mem = 2000
secpersite = 8 #max seconds per site - double what's expected for cf
time = secpersite*500
td = seconds_to_period(time)
time = sprintf('%02d:%02d:%02d', td@hour, minute(td), second(td))

l1 = "#!/bin/bash"
l4 = "#SBATCH --chdir=/scratch/ss3988/UVAFME"
l5 = "#SBATCH --mail-type=ALL"
l6 = paste0("#SBATCH --mem=",mem)
l7 = paste0("#SBATCH --time=",time)

f = "bsp"
for(m in mgmts) {
  for(deg in CC) {
    name = paste0("bsp_",m,"_",deg,"degCC")
    l3 = paste0("#SBATCH -o logs/site_mgmt/bsp/",name,"_err.log")
    l8 = paste0("name=file_lists/site_mgmt/bsp/file_list_",name,".txt")
    l9 = paste0("srun ",model," $name")
    
    lines = c(l1,l3,l4,l5,l6,l7,l8,l9)
    jobfile = paste0(jobdir, name, ".slurm")
    fileConn <- file(jobfile)
    writeLines(lines, fileConn)
    close(fileConn)
    joblist = c(joblist, paste0("sbatch ",jobfile))
  }
}

writeLines(joblist, "jobfiles/site_mgmt/bsp/joblist.txt")


