#!/usr/bin/env Rscript

#
# #########################
# Purpose: Create and move files for better batch scripting
# Author: Adrianna C. Foster
# Date: May, 2021
# R version 4.0.2 (2020-06-22) 'Taking Off Again'
# #########################
# Updated Jan, 2023 by Shelby Sundquist
# #########################
# #########################
# Input format: csv
# Output format: csv
# #########################

## Make sure we are searching the correct libraries
library(data.table)
library(dplyr)

#setwd("") # do this!
source('BatchScriptingFunctions.R')

## We have command-line arguments
args <- 'batch_params/batch_params_cf_r45.txt' #commandArgs(TRUE)
args <- c(args, 'batch_params/batch_params_cf_r85.txt') #commandArgs(TRUE)
args <- c(args, 'batch_params/batch_params_bau.txt')
args <- c(args, 'batch_params/batch_params_prj.txt')

wd = getwd() 

for(f in 1:4) {
## Parameter text file
param_file <- paste0(wd, "/",  as.character(args[f]))

## Read in parameters
paramValues <- get_params(param_file)

setwd(paramValues$mainDir)

## Check for errors
paramValues <- checkParams(paramValues)
units <- read.csv(paste0(paramValues$inputAll, paramValues$units))
paramValues$units = unique(units$unit)

## Create input and output folders
createDirs(paramValues)

## read in main sitelist
slist <- read.csv(paste0(paramValues$inputAll, paramValues$slist), stringsAsFactors = F)

if (!(paramValues$gcm == "")){
    gcm <- fread(paste0(paramValues$inputAll,paramValues$gcm), stringsAsFactors = F)
    gcm_lightning <- fread(paste0(paramValues$inputAll,paramValues$gcm_lightning), stringsAsFactors = F)
    batch.gcm = T
} else {
    batch.gcm = F
}

slist <- filter(slist, siteID %in% gcm$site & 
                    siteID %in% gcm_lightning$site)

## Subset sitelists and write to separate files/folders

sites = fread(paste0(paramValues$inputAll, paramValues$sites))

batchSitelists(slist = slist, inputMain = paramValues$inputMain, units = paramValues$units, sites = sites)

batchClimGCM(slist = slist, gcm = gcm, gcm_lightning = gcm_lightning, 
                units = paramValues$units, batch.gcm = batch.gcm, climDir = paramValues$climDir,
             sites = sites)

## Create file lists
  createFilelists(flistdir = paramValues$flistdir, inputAll = paramValues$inputAll, 
    inputMain =  paramValues$inputMain, outputMain =  paramValues$outputMain, units = paramValues$units, 
      rtdir = paramValues$rtdir, mgmtScenario = paramValues$mgmtscenario)

}

