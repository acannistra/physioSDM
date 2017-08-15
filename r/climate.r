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
    flog.info("\tOpening file %s", wc_file)
    wc_rast = raster(wc_file)
    occurrences[c(as.character(biovar))] = extract(wc_rast, occurrences[c(lonCol, latCol)])
  }
  
  return(occurrences)
}

getRasterStack_worldclim <- function(biovars, worldclimDir, extent=NULL, divide=NULL){
  brick = stack()
  for(biovar in biovars){
    wc_file = getWorldClimFilepath(biovar, worldclimDir)
    flog.info("\tOpening file %s", wc_file)
    wc_rast = raster(wc_file)
    if(!is.null(divide)){
      if(biovar %in% divide){
        wc_rast = wc_rast / 10
      }
    }
    names(wc_rast) = as.character(biovar)
    if(is.null(extent)){
      brick = stack(brick, wc_rast)
    } else {
      brick = stack(brick, crop(wc_rast, extent))
    }
  }
  return(brick)
}
####
#### WorldClim/BioCLIM accessors 
####

getWorldClimFilepath <- function(biovar, worldclimDir){
  if(!endsWith(worldclimDir, "/")) {
    stop(flog.error("worldclimDir needs to end with '/' character"))
    return("")
  }
  if(!is.numeric(biovar)){
    # if biovar is like "bio10" or "bio2"
    biovarid = as.integer(substring(biovar, 4))
  } else  {
    biovarid = biovar
  }
  filename = Sys.glob(sprintf("%s*%02d.tif", worldclimDir, biovarid))
  if(length(filename) == 0){
    stop(flog.error("No file for biovar %s in %s, aborting!", biovar, worldclimDir))
  }
  return(filename)
}

