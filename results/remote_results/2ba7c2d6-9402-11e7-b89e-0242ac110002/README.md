This experiment was run on all species, on an EC2 instance (c4.xlarge) with limited occurrence information (as a result of memory errors) to 3000 occurrences per species. See parameters.json for details on the experiment. 
// Parameter file for main.R program, defining 
// variables useful in experimentation. See documentation. 
// All paths must be absolute.
{
    "species" : "all", //list of specific species or "all"
    "max_occurrences": 3000,
    "chunksize" : 10,
    "time_periods" : ["1970:2009"], //past first, then future. 
    "bioclim_source" : "worldclim",
    "bioclim_dir" : "/projdir/data/climate/climate/current/",
    "biovars" : ["bio1", "bio3", "bio5", "bio6", "bio7", "bio12", "bio15"],    
    "prior_type" : "sigmoid",
    "phys_datafile" : "/projdir/data/physiology/Sundayetal_thermallimits_PNAS_2014.csv",
    "min_occurrence_threshold" : 100, //minimum number of GBIF occurrences necessary to build model.,
    "train_fraction" : 0.75,
    "plotting" : true, 
    "predict" : false,
    "loglevel" : "debug"
}
