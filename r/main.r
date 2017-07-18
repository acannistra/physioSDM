## Physiology in Species Distribution Models
## Tony Cannistra & Lauren Buckley, 2017
## ~ Main Experiment.

suppressPackageStartupMessages({
  library(futile.logger)
  library(jsonlite)
  library(dplyr)
  source('occurrences.R') ## for GBIFOccurrences and absence generation code.
  source("climate.R") ## for worldclim extraction code. 
})

PARAMETER_FILE = "../parameters.json" ## parameter file must contain absolute file paths. 

## 1) load parameters for experimentation 
PARAMS = fromJSON(PARAMETER_FILE)
flog.info("Loaded parameters from file '%s'", PARAMETER_FILE)

## 2) Build Species Experimentation Table
all_species = data.frame(species_name=PARAMS$species) %>%
                mutate(species_name = tolower(species_name))
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

testSpecies = function(speciesData){
  #theseOccs = GBIFOccurrences(speciesData[c('species_name')], currentTimePeriod_start, currentTimePeriod_end)
  flog.error("not implemented (species: %s)", speciesData[c('species_name')])
  ## get absence points
  ## connect climate data
  ## train model
  ## test model with holdout set
  ## get future climate data
  ## test model with future climate predictions -- work out details there. 
  return(NA)
}

apply(all_species, 1, testSpecies)
