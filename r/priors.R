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

gensigmoid <- function(x, low, high, rate, v, origin) {
  # [Generalized Sigmoid function.](https://en.wikipedia.org/wiki/Generalised_logistic_function)
  return(low + ((high-low)/(1+(rate*exp((x-origin))))^(1/v)))
}

buildPrior <- function(type, physioData, tminEnvCol='bio1', tmaxEnvCol='bio1'){
  ## assembles prior function given a type, a dataframe containing 'tmin' and 'tmax' for a single 
  ## species, and the corresponding covariate column names to which tmin and tmax priors are applied.
  ## returns a partially-applied function. 
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
        return(partial(sigmoid.range, 
                       tminEnvCol=tminEnvCol, tmaxEnvCol=tmaxEnvCol, 
                       tmin=tmin, tmax=tmax))
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
  result = ifelse(env<tmax, 0.5, gensigmoid(env, 0.1, 0.5, 5.5, 2.5, tmax))
  return(result)
}
sigmoid.tmin <- function(env, tmin, tminEnvCol) {
  env = env[,c(tminEnvCol)]
  result = ifelse(env>tmin, 0.5, gensigmoid(env, 0.5, 0.1, 5.5, 2.5, tmin))
  return(result)
}
sigmoid.range = function(env, tmax, tmin, tmaxEnvCol, tminEnvCol){
  result = c()
  
  evaluate_row = function(row){
    tmin_e_value = row[c(tminEnvCol)]
    tmax_e_value = row[c(tmaxEnvCol)]
    if (is.na(tmin_e_value) || is.na(tmax_e_value)){
      result = c(result, NA)
    } else if (tmin_e_value < tmin){
      result = c(result, gensigmoid(tmin_e_value, 0.5, 0.1, 5.5, 2.5, tmin))
    } else if (tmax_e_value > tmax){
      result = c(result, gensigmoid(tmax_e_value, 0.1, 0.5, 5.5, 2.5, tmax))
    } else{
      result = c(result, 0.5)
    }
  }
  apply(env, 1, evaluate_row)
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

  