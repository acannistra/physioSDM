suppressPackageStartupMessages({
  library(dismo)
  library(GRaF)
  library(dplyr)
  library(caret)
})

## Physiology in Species Distribution Models
## Tony Cannistra & Lauren Buckley, 2017
## SDM Helper Functions


#### Model Building Utilities

# create GRAF SDM model instance from dataframe of presence/absence occurrence
# information with covariates at each point. optional optimiziation (opt) and prior 
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

# Train MaxEnt model from dataframe of presence/absence occurrence
# information with covariates at each point
buildMaxEnt <-function(climPresAbs){
  labels = climPresAbs$presence
  covs   = dplyr::select(climPresAbs, -c(presence))
  return(maxent(covs, labels))
}


#### training and testing utilities

# create training/test data split.
trainTestSplit <- function(occurrences, trainFraction=0.75){
  split = list()
  indices = createDataPartition(1:nrow(occurrences), p=trainFraction)$Resample1
  split$train = occurrences[indices,]
  split$test  = occurrences[-indices,]
  return(split) 
}

#### Model Evaluation Functions

# compute AUC and ROC for a given GRAF model, covariates, and presence/absence labels.
graf_auc_roc = function(model, test_covs, test_labels){
  results = list()
  prob = data.frame(predict(model, test_covs))
  prob = prob$posterior.mode
  pred = prediction(prob, test_labels)
  auc  = performance(pred, measure='auc')
  auc = auc@y.values[[1]]
  results$auc = auc
  roc = performance(pred, measure='tpr', x.measure='fpr')
  results$roc = roc
  return(results)
}

# compute AUC and ROC for a given MaxEnt model, covariates, and presence/absence labels
maxent_auc_roc = function(model, test_covs, test_labels){
  results = list()
  prob = data.frame(predict(model, test_covs))
  pred = prediction(prob, test_labels)
  auc  = performance(pred, measure='auc')
  auc = auc@y.values[[1]]
  results$auc = auc
  roc = performance(pred, measure='tpr', x.measure='fpr')
  results$roc = roc
  return(results)
}