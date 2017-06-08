## prism_multiyear
# utlility functions for downloading multiple years' worth of PRISM climate data. 
library(prism)

prism_multiyear <- function(variable, years, dest = "~/prismtmp", months = 1:12){
	options(prism.path = dest)
	prism::get_prism_monthlys(variable, years, months, keepZip=F)
}