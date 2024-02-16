
#
# #########################
# Purpose: Create slurm jobs and txt file with lines to submit all jobs 
# Author: Shelby W. Sundquist
# Date: Jan, 2023
# R version 4.1.2
# #########################
# #########################
# Input format: csv
# Output format: txt, slurm 
# #########################

library(data.table)
library(dplyr)
library(lubridate)

#Specify working directory!
wd = ""

sites = fread("input_data/TVSF/clims_testing/main/UVAFME2018_site.csv")
mgmt = fread("input_data/TVSF/clims_testing/main/UVAFME2018_active_management.csv")

munits = unique(mgmt$unit)
nmunits = unique(sites$unit)
temp = as.data.frame(summary(as.factor(sites$unit)))
nunits = temp; nunits$unit = rownames(temp); colnames(nunits) = c("n","unit")


model = "UVAFME_init"
jobdir = "jobfiles/TVSF/clims_testing/"
flistdir = c("file_lists/TVSF/clims_testing/")

mscen = c("cf")
cscen = c("hist","gcm45","gcm85")
secpersite = 35*2 # estimate of max seconds to run each site. It shouldn't take longer than this and might take less than half the time for sites in small units
mempersite = 700 #0.7GB - sites should average 0.5GB each. For when running management ON

l1 = "#!/bin/bash"
l4 = paste0("#SBATCH --chdir=", wd)
l5 = "#SBATCH --mail-type=ALL"

#all units
scriptnames = c()

for(m in mscen[1]) {
  for(c in cscen) {
    for(u in nmunits) {
      name = paste0(m,"_",c,"_",u)
      time = secpersite*nunits[nunits$unit == u,]$n
      td = seconds_to_period(time)
      time = sprintf('%02d:%02d:%02d', td@hour, minute(td), second(td))
      mem = 2000 # no big reqs
      
      l2 = paste0("#SBATCH -o logs/TVSF/clims_testing/",name,"_out.log")
      l3 = paste0("#SBATCH -o logs/TVSF/clims_testing/",name,"_err.log")
      l6 = paste0("#SBATCH --mem=",mem)
      l7 = paste0("#SBATCH --time=",time)
      l8 = paste0("name=file_lists/TVSF/clims_testing/",m,"/",c,"/file_list_",u,".txt")
      l9 = paste0("srun ",model," $name")
      lines = c(l1,l2,l3,l4,l5,l6,l7,l8,l9)
      
      jobfile = paste0(jobdir, name, ".slurm")
      scriptnames = c(scriptnames, paste0("sbatch ", jobfile))

      fileConn <- file(jobfile)
      writeLines(lines, fileConn)
      close(fileConn)
    }
  }
}

fileConn <- file(paste0(jobdir, "slurmcalls.txt"))
writeLines(scriptnames, fileConn)
close(fileConn)

# Not sure why, but the cf_7C sites are all running out of memory. Raising these to 4000 by manually editing files
