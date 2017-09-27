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
  
  flog.info("building prior of type %s with tmin:%f (type: %s) and tmax:%f (type: %s)", type,
            tmin, physioData[c('tmin_metric')],
            tmax, physioData[c('tmax_metric')])
  flog.warn("only using critical thermal threshold data.")
  if(type == 'sigmoid'){
    if(!is.na(tmin) && !is.na(tmax)){
      ## we've got both tmin and tmax
      flog.info("\tBoth tmin and tmax present.")
      if(physioData[c('tmin_metric')] == 'crit' && physioData[c('tmax_metric')] == 'crit'){ ## both critical
        return(partial(sigmoid.neutralZone, 
                       tminEnvCol=tminEnvCol, tmaxEnvCol=tmaxEnvCol, 
                       neutralMin=tmin, neutralMax=tmax))
      } else if (physioData[c('tmin_metric')] != 'crit' &&physioData[c('tmax_metric')] == 'crit') {
        ## tmin is lethal, but tmax is critical, only tmax:
        flog.info("Only tmax used (tmin lethal).")
        return(partial(sigmoid.tmax,
                       tmaxEnvCol=tmaxEnvCol,
                       tmax=tmax))
      } else if(physioData[c('tmin_metric')] == 'crit' &&physioData[c('tmax_metric')] != 'crit'){
        # tmax is lethal, only tmin
        flog.info("Only tmin used (tmax lethal).")
        return(partial(sigmoid.tmin,
                       tminEnvCol=tminEnvCol,
                       tmin=tmin))
      }
    } else if (is.na(tmax) && (!is.na(tmin) && physioData[c('tmin_metric')] == 'crit')){
      flog.info("\tOnly tmin present.")
      ## we've got tmin but no tmax
      return(partial(sigmoid.tmin,
                     tminEnvCol=tminEnvCol,
                     tmin=tmin))
    } else if (is.na(tmin) && (!is.na(tmax) && physioData[c('tmax_metric')] == 'crit')){
      ## we've got tmax but no tmin
      flog.info("\tOnly tmax present.")
      return(partial(sigmoid.tmax,
                     tmaxEnvCol=tmaxEnvCol,
                     tmax=tmax))
    } else {
      flog.error("prior type %s not supported or only lethal physiology data available.", type)
      return(NULL)
    }
  } else if (type == 'thresh') {
    if(!is.na(tmin) && !is.na(tmax)){
      ## we've got both tmin and tmax
      flog.info("\tBoth tmin and tmax present.")
      if(physioData[c('tmin_metric')] == 'crit' && physioData[c('tmax_metric')] == 'crit'){ ## both critical
        return(partial(thresh.plateau, 
                       tminEnvCol=tminEnvCol, tmaxEnvCol=tmaxEnvCol, 
                       neutralMin=tmin, neutralMax=tmax))
      } else if (physioData[c('tmin_metric')] != 'crit' &&physioData[c('tmax_metric')] == 'crit') {
        ## tmin is lethal, but tmax is critical, only tmax:
        flog.info("Only tmax used (tmin lethal).")
        return(partial(thresh.tmax,
                       tmaxEnvCol=tmaxEnvCol,
                       tmax=tmax))
      } else if(physioData[c('tmin_metric')] == 'crit' &&physioData[c('tmax_metric')] != 'crit'){
        # tmax is lethal, only tmin
        flog.info("Only tmin used (tmax lethal).")
        return(partial(thresh.tmin,
                       tminEnvCol=tminEnvCol,
                       tmin=tmin))
      }
    } else if (is.na(tmax) || (!is.na(tmin) && physioData[c('tmin_metric')] == 'crit')){
      ## we've got tmin but no tmax
      flog.info("\tOnly tmin present.")
      return(partial(thresh.tmin,
                     tminEnvCol=tminEnvCol,
                     tmin=tmin))
    } else if (is.na(tmin) || (!is.na(tmax) && physioData[c('tmax_metric')] == 'crit')){
      ## we've got tmax but no tmin
      flog.info("\tOnly tmax present.")
  
      return(partial(thresh.tmax, 
                     tmaxEnvCol=tmaxEnvCol,
                     tmax=tmax))
    } else{
      stop(flog.error("prior type %s not supported or only lethal physiology data available.", type))
    }
  } else {
    stop(flog.error("prior type %s not supported or only lethal physiology data available.", type))
  }
}

####
#### Sigmoid functions with plateaus.
####

sigmoid.tmax <- function(env, tmax, tmaxEnvCol){
  env = env[,c(tmaxEnvCol)]
  result = ifelse(env<tmax, 0.5, exp(-(env-tmax)/5)-0.5)
  result[result <= 0] = 0.01
  result[result > 1] = 1
  return(result)
}
sigmoid.tmin <- function(env, tmin, tminEnvCol) {
  env = env[,c(tminEnvCol)]
  result = ifelse(env>tmin, 0.5, 0.5-exp(-(env-(tmin)/99000)))
  result[result <= 0] = 0.01
  result[result > 1] = 1
  return(result)
}
sigmoid.neutralZone <- function(env, tminEnvCol, tmaxEnvCol, neutralMin, neutralMax){
  tmin_prb = sigmoid.tmin(env, neutralMin, tminEnvCol)
  tmax_prb = sigmoid.tmax(env, neutralMax, tmaxEnvCol)
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

  