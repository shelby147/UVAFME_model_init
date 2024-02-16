#!/usr/bin/env Rscript

#
# #########################
# Purpose: Functions for batch scripting with UVAFME
# Author: Adrianna C. Foster
# Date: May, 2021
# R version 4.0.2 (2020-06-22) 'Taking Off Again'
# #########################
# Updated Jan, 2023 by Shelby Sundquist
# R Version 4.2.1 (2022-06-23) 'Funny-Looking Kid'
# #########################
# #########################
# Input format: csv
# Output format: csv
# #########################


createDirs <- function(paramValues, type = 'all') {
    
    ## Create input and output folders
    if (type == 'all') {
        createInputs(paramValues$inputMain, paramValues$units, "in")
        createInputs(paramValues$outputMain, paramValues$units, "out")
    } else if (type == 'input') {
        createInputs(paramValues$inputMain, paramValues$units, "in")
    } else if (type == 'output') {
        createInputs(paramValues$outputMain, paramValues$units ,"out")
    }
    
}

createInputs <- function(inputMain, units, type) {
    
    ## Create input folders
    if (!file.exists(inputMain)){
        dir.create(inputMain)
    }
    
    input_prefix <- 'unit_'

    for (i in units){
        
      if(type == "in") {
          inputDir <- paste0(inputMain, input_prefix, i)
          
          if (!file.exists(inputDir)){
              dir.create(inputDir)
          }
      } else if (type == "out") {
        outputDir <- paste0(inputMain, "hist/")
        if (!file.exists(outputDir)){
          dir.create(outputDir)
        }
        outputDir <- paste0(inputMain, "hist/", input_prefix, i)
        if (!file.exists(outputDir)){
          dir.create(outputDir)
        }
        
        outputDir <- paste0(inputMain, "gcm45/")
        if (!file.exists(outputDir)){
          dir.create(outputDir)
        }
        outputDir <- paste0(inputMain, "gcm45/", input_prefix, i)
        if (!file.exists(outputDir)){
          dir.create(outputDir)
        }
        
        outputDir <- paste0(inputMain, "gcm85/")
        if (!file.exists(outputDir)){
          dir.create(outputDir)
        }
        outputDir <- paste0(inputMain, "gcm85/", input_prefix, i)
        if (!file.exists(outputDir)){
          dir.create(outputDir)
        }
        
        
      } else {print("Not sure if we're creating inputs or outputs")}
        
    }
    
}

createOutputs <- function(outputMain, units) {
    
    ## Create input folders
    if (!file.exists(outputMain)){
        dir.create(outputMain)
    }
    
    output_prefix <- 'unit_'
    
    for (i in units){
        
        outputDir <- paste0(outputMain, output_prefix, i)
        
        if (!file.exists(outputDir)){
            dir.create(outputDir)
        }
        
    }
    
}

batchSitelists <- function(slist, inputMain, units, sites) {
  
    # Units correspond to TVSF management units 

    ## Take input sitelist and GCMs, and subset out to n files/folders
    
    # ## Number of sites per node
    # numSites <- nrow(slist)
    # nodesPerSite <- ceiling(numSites/nodes)
    # 
    # ## Add to slist file
    # nodelist <- rep(seq(1, nodes), each = nodesPerSite)
    # nodelist <- nodelist[1:numSites]
    # slist$node <- nodelist
  
    ## Create main input dir if we need to
    if (!file.exists(inputMain)){
        dir.create(inputMain)
    }
    
    ## Input dir prefix
    input_prefix <- 'unit_'
    
    for (i in units){
        
        ## Create input dir if we need to 
        inputDir <- paste0(inputMain, "/", input_prefix, i)
        if (!file.exists(inputDir)){
            dir.create(inputDir)
        }
        
        ## Subset site list
        temp = sites[sites$unit == i, ]$site
        sub = slist[slist$siteID %in% temp,]
        
        ## Write to file
        write.csv(sub, paste0(inputDir, '/UVAFME2018_sitelist.csv'), row.names = F, 
            quote = F, na = '')
        

        
    }
    
}

batchClimGCM <- function(slist, units, gcm = NA, gcm_lightning = NA, batch.gcm = F, 
                         climDir, sites) {
  
  if (!file.exists(climDir)){
    dir.create(climDir)
  }
  
  ## Input dir prefix
  input_prefix <- 'unit_'
  

  if (batch.gcm) {
    for (i in units){
      
          sub <- sites[sites$unit == i,]
      
          sub.gcm <- subset(gcm, gcm$site %in% sub$site)
          sub.gcm <- dplyr::arrange(sub.gcm, site, year)
          
          ## Create input dir if we need to 
          inputDir <- paste0(climDir, input_prefix, i)
          if (!file.exists(inputDir)){
            dir.create(inputDir)
          }
          
          ## Write to file
          write.csv(sub.gcm, paste0(inputDir, '/UVAFME2018_climate_GCM.csv'),
                    row.names = F, quote = F, na = '')

          sub.gcm.lightning <- subset(gcm_lightning, gcm_lightning$site %in% sub$site)
          sub.gcm.lightning <- dplyr::arrange(sub.gcm.lightning, site, year)

          ## Write to file
          write.csv(sub.gcm.lightning, paste0(inputDir, '/UVAFME2018_lightning_GCM.csv'),
              row.names = F, quote = F, na = '')
    }
  }
  
  
}

createFilelists <- function(flistdir, inputAll, inputMain, outputMain, units, mgmtScenario, rtdir){

    ## Create and write file lists
    
    ## Create filelist directory if it doesn't already exist
    if(!file.exists(flistdir)){
        dir.create(flistdir)
    }
  
    # Subdirectories for filelist by climate scenario
    flistdirhist = paste0(flistdir, "hist")
    if(!file.exists(flistdirhist)){
      dir.create(flistdirhist)
    }
    flistdirr45 = paste0(flistdir, "gcm45")
    if(!file.exists(flistdirr45)){
      dir.create(flistdirr45)
    }
    flistdirr85 = paste0(flistdir, "gcm85")
    if(!file.exists(flistdirr85)){
      dir.create(flistdirr85)
    }
  
  
    ## Input dir prefix
    input_prefix <- 'unit_'
    output_prefix <- 'unit_'
    
    inputAll <- (inputAll)
    
    # Create filelists
    
    # Historical climate
    for (i in units){
        
        ## Name of file list
        flist_file <- paste0(flistdirhist, '/file_list_', i, '.txt')
        
        inDir <- (paste0(inputMain,input_prefix, i))
        outDir <- (paste0(outputMain, "hist/", output_prefix,  i))
        rt <- paste0(rtdir, mgmtScenario, "/hist")
        
        ## Create lines
        l1 <- '&filenames'
        l2 <- paste0("input_directory='", inputAll, "'") # this is only rangelist
        l3 <- paste0("output_directory='", outDir, "'")
        l4 <- paste0("climate_directory='", inputAll, "'")
        l5 <- paste0("site_directory='", inputAll, "'")
        l6 <- paste0("sitelist_directory='", inDir, "'")
        l7 <- paste0("GCM_directory='", inputAll, "'") # doesn't matter for this one; no GCM being run
        l8 <- paste0("rt_directory='", rt, "'")
        l9 <- paste0("speclist_directory='", inputAll, "'")
        l10 <- '/'
        
        fileConn <- file(flist_file)
        flines <- c(l1, l2, l3, l4, l5, l6, l7, l8, l9, l10)
        writeLines(flines, fileConn)
        close(fileConn)
      
    }
    
    # GCM 4.5
    for (i in units){
      
      ## Name of file list
      flist_file <- paste0(flistdirr45, '/file_list_', i, '.txt')
      
      inDir <- (paste0(inputMain , input_prefix, i))
      outDir <- (paste0(outputMain, "gcm45/", output_prefix,  i))
      rt <- paste0(rtdir,mgmtScenario, "/gcm")
      GCMdir <- paste0("input_data/TVSF/clims_testing/","R45/",input_prefix ,i)
      
      ## Create lines
      l1 <- '&filenames'
      l3 <- paste0("output_directory='", outDir, "'")
      l6 <- paste0("sitelist_directory='", inDir, "'")
      l7 <- paste0("GCM_directory='", GCMdir, "'") 
      l8 <- paste0("rt_directory='", rt, "'")
      l10 <- '/'
      
      fileConn <- file(flist_file)
      flines <- c(l1, l2, l3, l4, l5, l6, l7, l8, l9, l10)
      writeLines(flines, fileConn)
      close(fileConn)
      
    }
    
    # GCM 8.5
    for (i in units){
      
      ## Name of file list
      flist_file <- paste0(flistdirr85, '/file_list_', i, '.txt')
      
      inDir <- (paste0(inputMain , input_prefix, i))
      outDir <- (paste0(outputMain, "gcm85/", output_prefix,  i))
      rt <- paste0(rtdir,mgmtScenario, "/gcm")
      GCMdir <- paste0("UVAFME/input_data/TVSF/clims_testing/","R85/",input_prefix ,i)
      
      
      ## Create lines
      l1 <- '&filenames'
      l3 <- paste0("output_directory='", outDir, "'")
      l6 <- paste0("sitelist_directory='", inDir, "'")
      l7 <- paste0("GCM_directory='", GCMdir, "'") 
      l8 <- paste0("rt_directory='", rt, "'")
      l10 <- '/'
      
      fileConn <- file(flist_file)
      flines <- c(l1, l2, l3, l4, l5, l6, l7, l8, l9, l10)
      writeLines(flines, fileConn)
      close(fileConn)
      
    }
    
}

checkParams <- function(paramValues) {
    
    ## Check parameter values
    
    cwd <- getwd()
    
    ## Create main directory if we have one and if it doesn't exist yet
    if (paramValues$mainDir == "" | is.na(paramValues$mainDir)){
        ## We don't have anything here
        paramValues$mainDir <- './'
        print(paste0('No main directory specified, using current working directory - ',
            cwd))
    } else {
        paramValues$mainDir <- paste0(paramValues$mainDir, '/')
        if (!file.exists(paramValues$mainDir)){
            dir.create(paramValues$mainDir)
        }
    }
    
    ## Input Main
    if (paramValues$inputMain == "" | is.na(paramValues$inputMain)){
        ## We don't have anything here
        paramValues$inputMain <- 'input_data/'
        print('No main input directory specified, using default of input_data')
    } else {
        paramValues$inputMain <- paste0(paramValues$mainDir, paramValues$inputMain, '/')
    }
    
    ## Output Main
    if (paramValues$outputMain == "" | is.na(paramValues$outputMain)){
        ## We don't have anything here
        paramValues$outputMain <- 'output_data/'
        print('No main output directory specified, using default of output_data')
    } else {
        paramValues$outputMain <- paste0(paramValues$mainDir, paramValues$outputMain, '/')
    }
    
    ## flistdir
    if (paramValues$flistdir == "" | is.na(paramValues$flistdir)){
        ## We don't have anything here
        paramValues$flistdir <- 'file_lists/'
        print('No filelist directory specified, using default of file_lists')
    } else {
        paramValues$flistdir <- paste0(paramValues$mainDir, paramValues$flistdir, '/')
    }
    
    ## rtdir
    if (paramValues$rtdir == "" | is.na(paramValues$rtdir)){
      ## We don't have anything here
      print('No filelist directory specified, using default of file_lists')
    } else {
      paramValues$rtdir <- paste0(paramValues$mainDir, paramValues$rtdir, '/')
    }
    
    ## climDir
    if (paramValues$climDir == "" | is.na(paramValues$climDir)){
      ## We don't have anything here
      print('No filelist directory specified, using default of file_lists')
    } else {
      paramValues$climDir <- paste0(paramValues$mainDir, paramValues$climDir, '/')
    }
    
    
    ## Input all
    if (paramValues$inputAll == "" | is.na(paramValues$inputAll)){
        ## We don't have anything here
        print('No input directory specified -- stopping.')
        stop()
    } else {
        paramValues$inputAll <- paste0(paramValues$mainDir, paramValues$inputAll, '/')
        if(!file.exists(paramValues$inputAll)){
            print(paste0('Directory ', paramValues$inputAll, 
                'does not exist -- stopping.'))
        }
    }
    
    ## sitelist
    if (paramValues$slist == "" | is.na(paramValues$slist)){
        ## We don't have anything here
        print('Sitelist file not specified -- stopping.')
        stop()
    } 
    
    return(paramValues)
    
}

get_params <- function(param_file){
    
## Read in parameter file
batch_params <- readLines(con = param_file)

## Empty list of parameters
paramValues <- list(mainDir = NA, inputMain = NA, outputMain = NA, 
                    inputAll = NA, units = NA, rtdir = NA, slist = NA, sites = NA, flistdir = NA,
                    gcm = NA, gcm_lightning = NA, climDir = NA, mgmtscenario = NA)

## Fill with values from input file
params <- names(paramValues)
for (i in 1:length(params)){
  ex <- paste0(params[i], "=")
  paramValues[[i]] <- gsub(ex, '', grep(ex, batch_params, value = T))
}

return(paramValues)
    
}

compileFiles <- function(outdirs, fileName){
    
    ## Gather output files from a batch run
    
    out.dat <- data.frame()
    
    for (i in 1:length(outdirs)){
        
        units <- gsub('unit_', '', basename(outdirs[i]))
        fname <- paste0(outdirs[i], '/', fileName)
        sub <- read.csv(fname, stringsAsFactors = F)
        sub$units <- units
        
        out.dat <- rbind(out.dat, sub)
        
    }
    
    return(out.dat)
    
}

get_status <- function(paramValues, outdirs) {
    
    ## Read in main sitelist 
    slist <- read.csv(paramValues$slist, stringsAsFactors = F)
    num_total <- nrow(slist)
    
    num_done <- 0
    for (i in 1:length(outdirs)){
        out_table <- read.table(paste0(outdirs[1], '/site_log.txt'))
        num_done <- num_done + length(unique(out_table$V3))
    }
    
    return(list(num_done = num_done, num_total = num_total))

    
}


