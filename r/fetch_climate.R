### fetch_monthly_climate 
library(raster)
library(prism)
library(optparse)

## parse user options
option_list = list(
  make_option(c("-p","--points"), type='character', help="CSV containing lon/lat points."),
  make_option("--sep", help='CSV separator character [optional, default comma]', default=','),
  make_option(c("-t","--time_resolution")),
  make_option(c("-m","--monthcol"), help="Column name with month values, if -t month"),
  make_option(c("-y","--yearcol"), help="Column name with year values."),
  make_option(c("-c", "--climate"), default = 'prism'),
  make_option("--latcol", help='Column name with latitude values.', default='decimallatitude'),
  make_option("--loncol", help='Column name with longitude values.', default='decimallongitude')
  
)
parser = OptionParser(option_list = option_list)
args = parse_args(parser)
##

## load points
if(!is.null(args$points)){
  # load file
  pointsdf = read.csv(args$points, sep=args$sep)
  # transform year into integer
  pointsdf[, args$yearcol] = as.integer(pointsdf[, args$yearcol])
  print(paste0("Found ",nrow(points)," points in file '",args$points,"'."))
}

pointsdf = SpatialPointsDataFrame(pointsdf[,c(args$loncol, args$latcol)], pointsdf)
## calculate required years
years = sort(na.omit(unique(pointsdf[[args$yearcol]])))
if (length(years) > 5) { print("Warning: number of requested years greater than 5.")}