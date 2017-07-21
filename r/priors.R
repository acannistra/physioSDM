suppressPackageStartupMessages({
  library(pryr)
})

## Physiology in Species Distribution Models
## Tony Cannistra & Lauren Buckley, 2017
## ~ Prior Functions for Physiological Data Incorporation

####
#### Notes: 
#### for the future, we should consider the difference between 
#### critical and lethal tmin and tmax values and adjust prior
#### assumptions accordingly. Current code does not do this. 


####
#### helper to build physiological prior based on 
#### known physiological data. 
####

buildPrior <- function(type, physioData, tminEnvCol='bio1', tmaxEnvCol='bio1'){
  tmin = as.numeric(physioData[c('tmin')])
  tmax = as.numeric(physioData[c('tmax')])
  
  #flog.info("building prior of type %s with tmin:%f and tmax:%f", type, tmin, tmax)
  if(type == 'sigmoid'){
    if(!is.na(tmin) && !is.na(tmax) && physioData[c('tmin_metric')] == 'crit' && physioData[c('tmax_metric')] == 'crit'){ ##only use crit
      ## we've got both tmin and tmax
      return(partial(sigmoid.neutralZone, 
                     tminEnvCol=tminEnvCol, tmaxEnvCol=tmaxEnvCol, 
                     neutralMin=tmin, neutralMax=tmax))
    } else if (!is.na(tmin) && is.na(tmax) && physioData[c('tmin_metric')] == 'crit'){
      ## we've got tmin but no tmax
      return(partial(sigmoid.tmin,
                     tminEnvCol=tminEnvCol,
                     tmin=tmin))
    } else if (is.na(tmin) && !is.na(tmax) && physioData[c('tmax_metric')] == 'crit'){
      ## we've got tmax but no tmin
      return(partial(sigmoid.tmax,
                     tmaxEnvCol=tmaxEnvCol,
                     tmax=tmax))
    }
  } else if (type == 'thresh') {
    if(!is.na(tmin) && !is.na(tmax)){ ##only use crit
      ## we've got both tmin and tmax
      return(partial(thresh.plateau,
                     tminEnvCol=tminEnvCol, tmaxEnvCol=tmaxEnvCol,
                     tmin=tmin, tmax=tmax))
    } else if (!is.na(tmin) && is.na(tmax) && physioData[c('tmin_metric')] == 'crit'){
      ## we've got tmin but no tmax
      return(partial(thresh.tmin,
                     tminEnvCol=tminEnvCol,
                     tmin=tmin))
    } else if (is.na(tmin) && !is.na(tmax) && physioData[c('tmax_metric')] == 'crit'){
      ## we've got tmax but no tmin
      return(partial(thresh.tmax, 
                     tmaxEnvCol=tmaxEnvCol,
                     tmax=tmax))
    }
  } else {
    stop(flog.error("prior type %s not supported or only lethal physiology data available.", type))
  }
}

####
#### Sigmoid functions with plateaus.
####

sigmoid.tmax <- function(env, tmax){
  return(ifelse(env<tmax, 0.9, exp(-(env-tmax)/5)))
}
sigmoid.tmin <- function(env, tmin) {
  result = ifelse(env<tmin, 0.1, 1-exp(-(env-(tmin)/9000)))
}
sigmoid.neutralZone <- function(env, tminEnvCol, tmaxEnvCol, neutralMin, neutralMax){
  tmin_prb = sigmoid.tmin(env[,c(tminEnvCol)], neutralMin)
  tmax_prb = sigmoid.tmax(env[,c(tmaxEnvCol)], neutralMax)
  return(tmin_prb*tmax_prb)
}

####
#### threshold functions with plateaus
####

thresh.tmax <- function(env, tmax) { 
  return(ifelse(env < tmax, 1.0, 0))
}
thresh.tmin <- function(env, tmin) {
  return(ifelse(env > tmin, 1.0, 0))
}
thresh.plateau <- function(env, tminEnvCol, tmaxEnvCol, tmin, tmax){
  return(thresh.tmax(env[,tmaxEnvCol], tmax) * thresh.tmin(env[,tminEnvCol], tmin))
}

  