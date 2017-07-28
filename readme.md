# physioSDM
> Code and experimentation for the inclusion of physiological data in species distribution models. 

## Experiment
We're trying to determine whether including physiological information as a prior in species distribution models can improve or change future projections. To do this we're building a physiologically-informed SDM for a list of species for which we have physiological data, and comparing the results to a non-physiological gaussian-process based SDM and MaxEnt. We compare the models in both their ability to recover current distributions and in how they project into the future. 


## Architecture. 
This is an experimentation suite built in R upon the following packages, all of which are dependencies:

* `dismo`
* `GrAF`
* `raster`
* `dplyr`, `pryr` and `purrr` (from the `tidy verse`)
* `ROCR`
* `uuid`
* `futile.logger`
* `jsonlite`
* `rgbif`
* `caret`

The code is architected such that it takes as input a .csv file containing various physiological parameters for a list of species, where each row is a species and each column represents some physiological quantity. The code also accepts a JSON file of parameters which determine many runtime behaviors. The same experimentation is performed 

All experimentation code is located in the `/r` directory. The following is an explanation of the files therein: 

*	`main.r`: contains the main experimentation loop and code to process input parameters.
*	`occurrences.R`: responsible for fetching occurrence information from GBIF. 
* 	`climate.R`: contains code for extracting climatological information from already-downloaded climate data which is stored locally. 
*  `sdm.r`: contains code for constructing all of the models used in experimentation. 
*  `priors.R`: contains code for constructing physiological priors for the gaussian process model. 

## Running Experiments
more later.
