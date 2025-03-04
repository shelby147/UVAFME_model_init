#
# #########################
# Purpose: Check output logs. Write files with jobfile names and error messages for any scenarios with errors. If all jobs ran properly, no new output files will be written
# Author: Shelby W. Sundquist
# Date: June, 2023
# R version 4.1.2
# #########################
# #########################
# Input format: log (txt)
# Output format: txt
# #########################
# 

library(data.table)
library(stringr)
library(sjmisc)

# info for reading files
logDir = "logs/KLC/"
jobDir = "jobfiles/KLC/"
cscen = c("hist","gcm45","gcm85")

sites = fread("input_data/KLC/main/UVAFME2018_site.csv")


goodlastline = "   ============================================================================="
joberrors = c()
errorcodes = c()

fname = paste0(logDir, "hist.err.log")
dat = readLines(fname)
nsites = as.integer(str_extract(dat[6], "[0-9]+"))
sitesrun = length( dat[grepl("Running for site", dat)])
lastline = dat[length(dat)]
if(lastline == goodlastline & nsites == sitesrun) {
  # scenario ran correctly
} else {
  jobname = paste0("sbatch ",jobDir, c, "_", i, ".slurm")
  joberrors = c(joberrors, jobname) 
  errorcodes = c(errorcodes, fname, lastline, paste0(sitesrun, " out of ", nsites, " run"))
}

for(c in cscen[2:3]) {
  for(i in 1:20) {
    fname = paste0(logDir, c, "_", i, ".err.log")
    dat = readLines(fname)
    nsites = as.integer(str_extract(dat[6], "[0-9]+"))
    sitesrun = length( dat[grepl("Running for site", dat)])
    lastline = dat[length(dat)]
    if(lastline == goodlastline & nsites == sitesrun) {
      # scenario ran correctly
    } else {
      jobname = paste0("sbatch ",jobDir, c, "_", i, ".slurm")
     joberrors = c(joberrors, jobname) 
     errorcodes = c(errorcodes, fname, lastline, paste0(sitesrun, " out of ", nsites, " run"))
    }
  }
}

if(length(joberrors) == 0) {
  joberrors = c(joberrors, "No job errors! Ready to compile output")
}
if(length(errorcodes) == 0) {
  errorcodes = c(errorcodes, "No job errors! Ready to compile output")
  print("No job errors! Ready to compile output")
}


badjobs = paste0(jobDir, "failedjobs.txt")
fileConn <- file(badjobs)
writeLines(joberrors, fileConn)
close(fileConn)

badjobs = paste0(jobDir, "failedjobserrors.txt")
fileConn <- file(badjobs)
writeLines(errorcodes, fileConn)
close(fileConn)

