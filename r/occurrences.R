suppressPackageStartupMessages({
  library(dismo)
  library(rgbif)
})

## Physiology in Species Distribution Models
## Tony Cannistra & Lauren Buckley, 2017
## ~ Occurrence helper functions. 

####
#### GBIF Data Extraction:
####
GBIFOccurrences = function(taxon, minyear=NA, maxyear=NA, limit=100000) {
  flog.info("searching GBIF for occurrences for species %s....", taxon)
  if (is.na(minyear) & is.na(maxyear)) {
    occs = occ_search(scientificName = taxon, hasCoordinate = T, limit = 100000)
  } 
  else {
    occs = occ_search(scientificName = taxon,
                      hasCoordinate = T, 
                      limit = limit,
                      eventDate = paste(minyear, maxyear, sep = ','))
  }
  if (is.null(occs$meta$count)){
    return(NULL)
  }
  return(occs$data)
}

#### 
#### Absence Generation
####

generateCircleAbsences <- function(occurrences, radius, count, latCol='decimalLatitude', lonCol='decimalLongitude'){
  circles = circles(occurrences[,c(lonCol, latCol)], d=radius, lonlat=T)
  pseudoabs = spsample(circles@polygons, count, type='random', iter=100)
  colnames(pseudoabs@coords) <- c(lonCol, latCol)
  return(as.data.frame(pseudoabs@coords))
}

generateRandomAbsences_bbox <- function(bbox, count){
  flog.error("generateRandomAbsences_bbox not implemented")
  return(data.frame())
}

generateRandomAbsences_poly <- function(poly, count){
  flog.error('generateRandomAbsences_poly not implemented')
  return(data.frame())  
}


####
#### Utilities
####

plotOccurrences = function(occs){
  library(ggplot2)
  library(ggmap)
  lon_range = extendrange(occs$decimalLongitude)
  lat_range = extendrange(occs$decimalLatitude)
  zoom = calc_zoom(lonrange, latrange)
  
  lon_center = mean(occs$decimalLongitude)
  lat_center = mean(occs$decimalLatitude)
  map = get_map(location = c(lon_center, lat_center), zoom)
  
  ggmap(map) + geom_point(
    aes(x = decimalLongitude, y = decimalLatitude), 
    data=occs, 
    color = 'blue', 
    size=1
  )
  
}

mergePresAbs <- function(pres, abs){
  if (ncol(pres) != ncol(abs)){
    flog.error("Presence and Absence dataframes do not have same number of columns.
               Be sure that they have followed the same processing pipeline (i.e if climate data was added
               to one set, be sure it is added to the other set before merging.")
    return()
  }
  pres[c('presence')] = 1
  abs[c('presence')]  = 0
  return(rbind(pres, abs))
}



