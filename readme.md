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

## Installation.

1. Clone this repository to your local computer. 
2. Acquire the datasets for your experimentation. More information on how to get our specific data will come with publication. 
4. Edit `parameters.json` file to reflect these data and your chosen parameters. 
3. Install dependencies. 
	1. `apt-get install libgeos-dev libgdal-dev default-jdk`. *Note that depending on your system you may need to install additional libraries. You'll know this after the next step; R will tell you what's missing.*
	2. Run `install.packages(c('dismo', 'GRaF', 'raster', 'tidyverse', 'ROCR', 'uuid', 'futile.logger', 'jsonlite', 'rgbif', 'caret', 'rgdal', 'rJava'))` in R. 
4. Attempt dry run by changing into `r/` directory and running `Rscript main.r`. You will likely need to download `maxent.jar` and install it. The program will instruct you on how to do this. 

 

## Running Experiments
more later.
