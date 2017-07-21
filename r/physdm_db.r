library(dismo)
library(raster)
library(dplyr)
library(GRaF)
library(RPostgreSQL)

.dbuser = 'postgres'
.dbpass = 'defaultpassword'
.dbhost = ''
.dbname = 'climate'
.dbtable  = 'bioclim_yearly'

dbConnection = dbConnect(dbDriver("PostgreSQL"),
                         db=.dbname,
                         host=.dbloc,
                         pass=.dbpass, 
                         user=.dbuser)


