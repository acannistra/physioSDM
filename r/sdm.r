suppressPackageStartupMessages({
  library(dismo)
  library(GRaF)
  library(dplyr)
  library(caret)
})

## Physiology in Species Distribution Models
## Tony Cannistra & Lauren Buckley, 2017
## SDM Helper Functions


####
#### Model Building Utilities
####

buildSDM <- function(climPresAbs, opt=F, prior=NULL){
  labels = climPresAbs$presence
  covs   = dplyr::select(climPresAbs, -c(presence))
  prior_name = as.character(quote(prior))
  if(is.null(prior)){
    if(!opt){
      flog.info("Training gaussian process model without prior and without optimization. ")
    } else{
      flog.info("Training and optimizing gaussian process model without prior.")
    } 
  } else {
      if(!opt){
        flog.info("Training gaussian process model with prior function '%s' without optimization.", prior_name)
      } else {
        flog.info("Training and optimizing gaussian process model with prior function '%s'", prior_name)
      }
  }
  model = graf(as.numeric(labels), covs, opt.l=opt, prior=prior)
  flog.info("...done.")
  return(model)
}

####
#### training and testing utilities
####

trainTestSplit <- function(occurrences, trainFraction=0.75){
  split = list()
  indices = createDataPartition(1:nrow(occurrences), p=trainFraction)$Resample1
  split$train = occurrences[indices,]
  split$test  = occurrences[-indices,]
  return(split) 
}
