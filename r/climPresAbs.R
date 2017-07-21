library(RPostgreSQL)
library(dict)
library(lubridate)
library(raster)
library(dplyr)
library(dismo)

#
# The purpose of this file is to define the function `climPresAbs`
# which exists to provide presence/absence data for a given database
# of survey occurrence information. 
#
# IMPORTANT: input occurrence information should be for a *SINGLE* taxon.
# All records will be considered to be from the same survey (i.e. no further
# data selection is completed herein aside from the omission of rows with
# missing data.)
#
# This function relies on a functional database containing 
# bioclim variables keyed on (lat, lon, year).
#


# Default climate database connection parameters
.dbName = "climate"
.dbLoc  = 'localhost'
.dbUser = 'tony'
.dbPass = ""
.tablename = "bioclim_geom"



db.climPresAbs <- function (occs, biovars, pseudoabsFrac = 1,
                            tablename = .tablename, dbname = .dbName, dbuser = .dbUser, dbpass = .dbPass, dbloc = .dbLoc){
  
  #
  # Given occurrence information, generates presence and absence
  # records with accompanying climate data (assumes local climate DB).
  #
  # returns a dataframe containing a "presence" column, lat, lon, 
  # and all bioclim vars.
  #
  # IMPORTANT: input occurrence information should be for a *SINGLE* taxon.
  # All records will be considered to be from the same survey (i.e. no further
  # data selection is completed herein aside from the omission of rows with
  # missing data.)
  dbConnection = dbConnect(dbDriver("PostgreSQL"),
                           db=dbname,
                           host=dbloc,
                           pass=dbpass, 
                           user=dbuser)
  
          
  num_pseudoabs = nrow(occs) * pseudoabsFrac
  bufs = circles(occs[,c("decimallongitude", "decimallatitude")], d=50000, lonlat=TRUE)
  abs  = spsample(bufs@polygons, num_pseudoabs, type='random', iter=100)

  ## get climate data
  ### must check if occurences spread multiple years
  occ_years = distinctYears(occs)
  climRasters = fetchDbClim.many_years(dbConnection, tablename, biovars, occ_years)
  
  ## assign data to each presence ...
  clim_Pres = occs %>%
              rowwise() %>%
              do(as.data.frame(extract(chooseRaster(climRasters, .$year), t(as.matrix(c(.$decimalLongitude, .$decimalLatitude)))))) %>%
              select(-year) ## remove "year" column.

  if (all(is.na(clim_Pres))) {
    return(NULL)
  }
  
  clim_Pres = data.frame(lon=occs$decimallongitude,
                         lat=occs$decimallatitude,
                         clim_Pres)

  ## ..and absence point. gets climate data from most frequently-seen year in occurrences. not sure if this is good.
  clim_Abs = as.data.frame(extract(chooseRaster(climRasters, mostFrequentYear(occs)), abs)) %>% select(-year)

        
  
  clim_Abs  = data.frame(lon=abs@coords[,'x'], lat=abs@coords[,'y'], clim_Abs)

  
  ## assign binary presence (1) absence (0) flags.
  presence = rep(1,dim(clim_Pres)[1])
  presence_temp = data.frame(presence, clim_Pres)
  presence = rep(0, dim(clim_Abs)[1])
  absence_temp = data.frame(presence, clim_Abs)
  
  
  
  ## and combine them. 
  clim_PresAbs = rbind(presence_temp, absence_temp)
  return(clim_PresAbs)
}

distinctYears <- function(df){
  #
  # for a given data frame with a "year" column, 
  # computes a list of distinct years represented
  #
  return(as.list(df %>% 
                 select(year) %>% 
                 na.omit() %>% 
                 distinct())$year)
}

mostFrequentYear <- function(df){
  #
  # for a given data frame with a "year column,
  # computes the year with the most number of 
  # corresponding records. 
  #
  return(as.list(df %>%
         group_by(year) %>% 
         count(sort = T) %>%
         top_n(1))$year)
}

chooseRaster <- function(rasters, year){
  #
  # given a dict() of rasters keyed on numeric years,
  # returns the raster corresponding to that year.
  # Key must exist, or else returns NULL.
  #
  return(rasters$get(year))
}

fetchDbClim.by_year <- function(connection, tablename, biovars, year){
  #
  # Retrieves biovars from db at connection corresponding to year
  # Returns a Raster.
  #
  
  queryFmt = "SELECT x, y, year, %s FROM %s WHERE year=%d;"
  query    = sprintf(queryFmt, paste(biovars, collapse=','), tablename, year)
  print(query)
  result   = data.frame(dbGetQuery(connection, query))
  if(nrow(result) == 0){
    stop("No results from database.")
  }
  return(rasterFromXYZ(result))
}

fetchDbClim.many_years <- function(connection, tablename, biovars, years){
  #
  # Retrieves biovars from db at connection corresponding to all 
  # years in year. 
  #
  # Returns a dictionary of rasters, keyed on year.
  #
  clim <- dict()
  for(year in years){
    clim[[year]] = fetchDbClim.by_year(connection, tablename, biovars, year)
  }
  
  return(clim)
}
