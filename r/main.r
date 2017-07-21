## Physiology in Species Distribution Models
## Tony Cannistra & Lauren Buckley, 2017
## ~ Main Experiment.

suppressPackageStartupMessages({
  library(futile.logger)
  library(jsonlite)
  library(dplyr)
  source('occurrences.R') ## for GBIFOccurrences and absence generation code.
  source("climate.R") ## for worldclim extraction code. 
  source("sdm.R") ## for gaussian process SDM code + helper functions.
  source("priors.R") ## for physiological priors. 
})

PARAMETER_FILE = "../parameters.json" ## parameter file must contain absolute file paths. 

## 1) load parameters for experimentation 
PARAMS = fromJSON(PARAMETER_FILE)
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


speciesExperiment = function(speciesData){
  flog.info("Starting experimentation with species %s", speciesData[c('species_name')])
  #thesePres = GBIFOccurrences(speciesData[c('species_name')], currentTimePeriod_start, currentTimePeriod_end)
  thesePres = GBIFOccurrences(speciesData[c('species_name')], currentTimePeriod_start, currentTimePeriod_end, limit = 100)
  if(is.null(thesePres)){
    flog.warn("No occurrences from GBIF for this species! Aborting experiment and continuing.")
    flog.info("*********species %s complete**********", speciesData[c('species_name')])
    
    return(NULL)
  }
  ## set up plotting parameters:
  plotrows = ceiling(length(PARAMS$biovars) / 3)
  plotcols = 3

  ## get absence points
  flog.info("Generating absence points using circle method.")
  theseAbs  = generateCircleAbsences(thesePres, 50000, nrow(thesePres))
  
  ## Merge presence and absence data. 
  all_points = mergePresAbs(thesePres[c('decimalLongitude', 'decimalLatitude')], theseAbs)
  flog.info("Merged %d points.", nrow(all_points))
  
  ## connect climate data
  flog.info("Assigning climate data using source %s located at %s.", PARAMS$bioclim_source, PARAMS$bioclim_dir)
  all_with_clim = assignPointData_worldclim(all_points, PARAMS$biovars, PARAMS$bioclim_dir)
  ### drop lat/lon and remove NAs. 
  all_with_clim = dplyr::select(all_with_clim, -c(decimalLongitude, decimalLatitude))
  completecases = complete.cases(all_with_clim)
  all_with_clim = data.frame(all_with_clim[completecases,])
  
  
  ## get training/test set (maybe)
  flog.warn("training/test separation NOT IMPLEMENTED")
  
  ## train model without prior
  model_noprior = buildSDM(all_with_clim, opt=T)
  par(mfrow=c(plotrows, plotcols))
  if(PARAMS$plotting) plot(model_noprior, data=F, jitter=0) ## TODO: There's something wrong with plotting, rugs dont work

  
  ## develop prior
  tmin = as.numeric(speciesData[c('tmin')])
  tmax = as.numeric(speciesData[c('tmax')])
  print(typeof(tmin))
  prior = buildPrior('sigmoid', tminEnvCol = 'bio1', tmaxEnvCol = 'bio1', speciesData)

  
  
  ## train model with prior
  
  
  model_withprior = buildSDM(all_with_clim, prior = prior, opt=T)
  plot.new()
  par(mfrow=c(plotrows, plotcols))
  if(PARAMS$plotting) plot(model_withprior, data=T, prior=T)
  
  ## test both models for current distirbution recovery 
  ## get future climate data
  ## test model with future climate predictions -- work out details there. 
  flog.info("*********species %s complete**********", speciesData[c('species_name')])
  return("not implemented")
}

apply(all_species, 1, speciesExperiment)
