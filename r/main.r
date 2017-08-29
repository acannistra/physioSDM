## Physiology in Species Distribution Models
## Tony Cannistra & Lauren Buckley, 2017
## ~ Main Experiment.

suppressPackageStartupMessages({
  library(futile.logger)
  library(jsonlite)
  library(dplyr)
  library(ROCR)
  library(purrr)
  library(uuid)
  source('occurrences.R') ## for GBIFOccurrences and absence generation code.
  source("climate.r") ## for worldclim extraction code. 
  source("sdm.r") ## for gaussian process SDM code + helper functions.
  source("priors.R") ## for physiological priors. 
})
PARAMETER_FILE = "../parameters.json" ## parameter file must contain absolute file paths. 
WORLDCLIM_SCALED_BIOCLIM_VARS = c("bio1", 
                                  "bio2", 
                                  "bio5", 
                                  "bio6", 
                                  "bio7", 
                                  "bio8", 
                                  "bio9", 
                                  "bio10", 
                                  "bio11")
OLD.PAR <- par(mar = c(0, 0, 0, 0))


## 1) load parameters for experimentation 
PARAMS = fromJSON(PARAMETER_FILE)
if (PARAMS$loglevel == 'debug'){
  flog.threshold(DEBUG)
}
SESSION_UUID = UUIDgenerate(use.time=T)
CODE_VERSION = system("git rev-parse --short HEAD", intern = TRUE)
flog.info("Session UUID: %s", SESSION_UUID)
flog.info("Code Version: %s", CODE_VERSION)
flog.info("Loaded parameters from file '%s'", PARAMETER_FILE)


## 2) Build Species Experimentation Table
if(PARAMS$species == 'all'){
  ## use all species with physiological information.
  physiology_temp = read.csv(PARAMS$phys_datafile)
  physiology_temp = physiology_temp %>% mutate(species_clean = tolower(gsub('_', ' ', species)))
  all_species = data.frame(species_name = physiology_temp$species_clean)
  rm(physiology_temp)
} else{ 
  ## use species listed in PARAMS$species
  all_species = data.frame(species_name=PARAMS$species) %>%
    mutate(species_name = tolower(species_name))
}
flog.info("Found %d species to analyze in parameter file.", nrow(all_species))



## 3) Add species physiological information (if available) to table.
flog.info("Reading physiology data from %s", PARAMS$phys_datafile)
physiology_data = read.csv(PARAMS$phys_datafile)

# transform names for consistency
physiology_data = physiology_data %>%
                    mutate(clean_name = tolower(gsub('_', ' ', species)))
# join species names with physiology information. 
all_species = all_species %>% 
              left_join(physiology_data, c("species_name" = "clean_name"))

numWithoutPhysio = length(which(apply(subset(all_species, select = -c(species, species_name)),
                                   1, function(x)all(is.na(x)))))

flog.info("Loaded physiology data for %d species (%d species with no data.)", nrow(all_species) - numWithoutPhysio, numWithoutPhysio)

## 4) Begin Experimentation
## 4a) define function for each species, to operate on each row of the all_species dataframe. 
currentTimePeriod_start = unlist(strsplit(PARAMS$time_periods[1], ":"))[1]
currentTimePeriod_end   = unlist(strsplit(PARAMS$time_periods[1], ":"))[2]

.FAILFILE = function(species, reason, file){
  write(sprintf("%s\t%s", species, reason), append=T, file=file)
}

speciesExperiment = function(speciesData){
  flog.info("Starting experimentation with species %s", speciesData[c('species_name')])
  #thesePres = GBIFOccurrences(speciesData[c('species_name')], currentTimePeriod_start, currentTimePeriod_end)
  thesePres = GBIFOccurrences(speciesData[c('species_name')], currentTimePeriod_start, currentTimePeriod_end, limit = 10000)
  if(is.null(thesePres)){
    flog.warn("No occurrences from GBIF for this species! Aborting experiment and continuing.")
    flog.info("*********species %s complete**********", speciesData[c('species_name')])
    .FAILFILE(speciesData[c('species_name')], "no occurrences.", failfile)
    return(NULL)
  } else if (nrow(thesePres) < PARAMS$min_occurrence_threshold){
    flog.warn("Number of GBIF occurrences (%d) less than min_occurrence_threshold parameter (%d)", nrow(thesePres), PARAMS$min_occurrence_threshold)
    flog.info("*********species %s complete**********", speciesData[c('species_name')])
    .FAILFILE(speciesData[c('species_name')], "minimum occurrence threshold not met", failfile)
    return(NULL)
  } else {
    flog.info("found %d occurrences.", nrow(thesePres))
  }

  ## set up plotting parameters for graf
  plotrows = ceiling(length(PARAMS$biovars) / 3)
  plotcols = 3
  plotwidth = 800
  plotheight = 600
  if(PARAMS$plotting){
    dir.create(sprintf("%s/plots/%s", results_dir, gsub(" ", "-", speciesData[c('species_name')])))
  }


  ## get absence points
  all_points = tryCatch({
    flog.info("Generating absence points using circle method.")
    theseAbs  = generateCircleAbsences(thesePres, 50000, nrow(thesePres))
  
    ## Merge presence and absence data. 
    all_points = mergePresAbs(thesePres[c('decimalLongitude', 'decimalLatitude')], theseAbs)
    flog.info("Merged %d points.", nrow(all_points))
    all_points
    }, error=function(e){
      flog.error(sprintf("failure in circles: %s", str(e$message)))
      return(NULL)
    })
  
  if(is.null(all_points)){
    .FAILFILE(speciesData[c('species_name')], "error in absence generation or merging", failfile)
    return(NULL)
  }
  ## connect climate data
  flog.info("Assigning climate data using source %s located at %s.", PARAMS$bioclim_source, PARAMS$bioclim_dir)
  all_with_clim = assignPointData_worldclim(all_points, PARAMS$biovars, PARAMS$bioclim_dir)
  ### drop lat/lon and remove NAs. 
  all_with_clim = dplyr::select(all_with_clim, -c(decimalLongitude, decimalLatitude))
  completecases = complete.cases(all_with_clim)
  all_with_clim = data.frame(all_with_clim[completecases,])
  
  
  ## get training/test set (maybe)
  split = trainTestSplit(all_with_clim, trainFraction=PARAMS$train_fraction)
  flog.info("Witholding %.2f%% of the data for training (%d occurrences)", 100*PARAMS$train_fraction, nrow(split$train))
  
  ## train model without prior
  model_noprior = tryCatch({
    buildSDM(split$train, opt=F)
  }, error = function(e){
    flog.error(str(e$message))
    return(NULL)
  })
  if (is.null(model_noprior)){
    .FAILFILE(speciesData[c('species_name')], "model without prior failed to build.", failfile)
    return(NULL)
  } 
  
  par(mfrow=c(plotrows, plotcols))
  if(PARAMS$plotting) {
    png(sprintf("%s/plots/%s/graf_responses.png", results_dir, gsub(" ", "-", speciesData[c('species_name')])),
        width=plotwidth,
        height=plotheight)
    par(mfrow=c(plotrows, plotcols))
    plot(model_noprior, data=F, jitter=0) ## TODO: There's something wrong with plotting, rugs dont work
    title(sub = sprintf("(code version: %s)", CODE_VERSION),
          cex.main = 2,   font.main= 2, col.main= "blue",
          cex.sub = 0.75, font.sub = 3, col.sub = "red")
    dev.off()
  }
  ## develop prior
  tmin = as.numeric(speciesData[c('tmin')])
  tmax = as.numeric(speciesData[c('tmax')])
  prior = tryCatch({
    buildPrior(PARAMS$prior_type, tminEnvCol = 'bio6', tmaxEnvCol = 'bio5', speciesData)
  }, error = function(e){
    return(NULL)
  })
  ## train model with prior
  model_withprior = NULL
  if(!is.null(prior)){
    model_withprior = tryCatch({
      buildSDM(split$train, prior = prior, opt=F)
    }, error = function(e){
      flog.error(str(e$message))
      return(NULL)
    })
  }
  
  if (is.null(model_withprior)){
    .FAILFILE(speciesData[c('species_name')], "model with prior failed to build", failfile)
    return(NULL)
  } 
  
  plot.new()
  if(PARAMS$plotting) {
    png(sprintf("%s/plots/%s/graf2_prior_responses.png", results_dir, gsub(" ", "-", speciesData[c('species_name')])), 
        width=plotwidth,
        height=plotheight)
    par(mfrow=c(plotrows, plotcols))
    plot(model_withprior, data=T, prior=T)
    title(sub = sprintf("(code version: %s)", CODE_VERSION),
          cex.main = 2,   font.main= 2, col.main= "blue",
          cex.sub = 1, font.sub = 3, col.sub = "red")
    dev.off()
  }
  
  ## test both models for current distirbution recovery 
  png(sprintf("%s/plots/%s/ROC.png", results_dir, gsub(" ", "-", speciesData[c('species_name')])),
      width=plotwidth,
      height=plotheight)
  par(mfrow=c(1,1))
  
  testCovs   = subset(split$test, select = -presence)
  testLabels = subset(split$test, select = presence)
  
  noprior_aucroc = graf_auc_roc(model_noprior, testCovs, testLabels)
  flog.info("No Prior ROC (AUC: %2f)", noprior_aucroc$auc)
  plot(noprior_aucroc$roc, col=rainbow(10), label='no prior')
  abline(a=0, b=1)

  prior_aucroc = graf_auc_roc(model_withprior, testCovs, testLabels)
  plot(prior_aucroc$roc, col='blue', add=T)

    
  flog.info("With Prior ROC (AUC: %2f)", prior_aucroc$auc)

  ## Compare against MaxEnt
  try(dyn.load("/Library/Java/JavaVirtualMachines/jdk1.8.0_111.jdk/Contents/Home/jre/lib/server/libjvm.dylib")) ##only on mac
  model_maxent = buildMaxEnt(split$train)
  maxent_aucroc = maxent_auc_roc(model_maxent, testCovs, testLabels)
  plot(maxent_aucroc$roc, col='green', add=T)

  legend(0.5, 0.2, c(sprintf("No Prior (AUC: %2f)", noprior_aucroc$auc),
                      sprintf("With Prior (AUC: %2f)", prior_aucroc$auc),
                      sprintf("MaxEnt (AUC: %2f", maxent_aucroc$auc)), 
         lty=c(1,1), 
         lwd=c(2.5,2.5),col=c("red","blue", "green"))
  title(main=sprintf("%s Model Performance",  speciesData[c('species_name')]), 
        sub = sprintf("(code version: %s)", CODE_VERSION),
        cex.sub = 1, font.sub = 3, col.sub = "red")
  dev.off()
  
  ## test model with future climate predictions -- work out details there. 
  ## get future climate data
  if(PARAMS$predict){
    testExtent = eval(parse(text=PARAMS$test_extent))
    futureClim = getRasterStack_worldclim(PARAMS$biovars, 
                                          PARAMS$future_bioclim_dir, 
                                          divide=union(PARAMS$biovars, WORLDCLIM_SCALED_BIOCLIM_VARS))
    points = data.frame(rasterToPoints(crop(futureClim, testExtent)))
    latLonFuture = subset(points, select=c(x, y))
    varsFuture = subset(points, select=-c(x, y))
    ans = cbind(latLonFuture, data.frame(predict(model_noprior, varsFuture))$posterior.mode)
    par(OLD.PAR)
    plot.new()
    plot(rasterFromXYZ(ans))
  }

  
  ## Complete. Return Results
  
  results = data.frame(species=speciesData[c('species_name')],
                       graf_cur_auc = noprior_aucroc$auc, 
                       phys_cur_auc = prior_aucroc$auc, 
                       maxent_cur_auc = maxent_aucroc$auc)
  
  flog.info("*********species %s complete**********", speciesData[c('species_name')])
  return(results)
}

### Run Experiment on All Species.

#### Create Results + Failures Storage Directories
results_dir = sprintf("../results/%s", SESSION_UUID)
dir.create(results_dir)
if(PARAMS$plotting){ dir.create(sprintf("%s/plots", results_dir))}

if(!is.null(PARAMS$chunksize)){
  chunk_id = 0
  chunks <- ggplot2::cut_interval(1:nrow(all_species), length=PARAMS$chunksize, labels=FALSE)
  for (chunk in unique(chunks)){
    these_species = all_species[which(chunks==chunk),]
    flog.info("Starting chunk %d (length: %d)", chunk_id, nrow(these_species))
    
    failures_fname = sprintf("%s/failures-chunk%d.txt", results_dir, chunk_id)
    results_fname = sprintf("%s/results-chunk%d.csv", results_dir, chunk_id)
    
    #### Copy parameters into runtime directory
    file.copy(PARAMETER_FILE, sprintf("%s/parameters.json", results_dir))
    
    #### Initialize Failure File
    failfile = file(failures_fname, 'w')
    
    #### RUN EXPERIMENT + Save results
    results = reduce(apply(these_species, 1, speciesExperiment), rbind)
    write.table(results, file=results_fname, sep=',', row.names=FALSE)
    chunk_id = chunk_id + 1
  }
}

failures_fname = sprintf("%s/failures.txt", results_dir)
results_fname = sprintf("%s/results.csv", results_dir)
dir.create(results_dir)
if(PARAMS$plotting){ dir.create(sprintf("%s/plots", results_dir))}

#### Copy parameters into runtime directory
file.copy(PARAMETER_FILE, sprintf("%s/parameters.json", results_dir))

#### Initialize Failure File
failfile = file(failures_fname, 'w')

#### RUN EXPERIMENT + Save results
results = reduce(apply(all_species, 1, speciesExperiment), rbind)
write.table(results, file=results_fname, sep=',', row.names=FALSE)
  