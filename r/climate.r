suppressPackageStartupMessages({
  library(raster)
})

## Physiology in Species Distribution Models
## Tony Cannistra & Lauren Buckley, 2017
## ~ Climate Data Helper Functions


####
#### Convenience wrappers
####

## extracts biovars for each occurrence from files in worldclimDir (which must end in /).

assignPointData_worldclim <- function(occurrences, biovars, worldclimDir, latCol='decimalLatitude',  lonCol = 'decimalLongitude'){
  for (biovar in biovars){
    wc_file = getWorldClimFilepath(biovar, worldclimDir)  
    wc_rast = raster(wc_file)
    occurrences[c(biovar)] = extract(wc_rast, occurrences[c(lonCol, latCol)])
  }
  return(occurrences)
}

####
#### WorldClim/BioCLIM accessors 
####

getWorldClimFilepath <- function(biovar, worldclimDir){
  if(!endsWith(worldclimDir, "/")) {
    flog.error("worldclimDir needs to end with '/' character")
    return("")
  }
  if(startsWith(biovar, "bio")){
    # if biovar is like "bio10" or "bio2"
    biovarid = as.integer(substring(biovar, 4))
  } else if (is.numeric(biovar)) {
    biovarid = biovar
  } else {
    flog.error("biovar variable is not numeric.")
  }
  file = Sys.glob(sprintf("%s*%02d.tif", worldclimDir, biovarid))
  print(file)
}
