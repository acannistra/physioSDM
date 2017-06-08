library(prism)
source("prism_multiyear.R")
## fetch_prism.R
## usage: fetch_prism.R {variables: comma separated} {timetype (day, month, normal)} {yearrange: YYYY:YYYY}

args = commandArgs(trailingOnly = T)
if (length(args) < 3){
	stop("usage: fetch_prism.R {variables: comma separated} {timetype (day, month, normal)} {yearrange: YYYY:YYYY}")
}
VARIABLES = strsplit(args[1], ",")[[1]]
TIMETYPE  = args[2]
START_Y   = strsplit(args[3], ":")[[1]][1]
END_Y	  = strsplit(args[3], ":")[[1]][2]
if(length(args) == 4){
	PRISMDIR = args[4]
} else {
	PRISMDIR = "~/prismtmp"
}

if (TIMETYPE == 'day'){
	stop("Day Not Implemented")
} else if (TIMETYPE == 'month') {
	for (var in VARIABLES){
		prism_multiyear(var, START_Y:END_Y)
	}
} else if (TIMETYPE == 'normal') {
	stop("Normals Not Implemented")
} else {
	stop(paste0("Time type ", TIMETYPE, " not supported."))
}






