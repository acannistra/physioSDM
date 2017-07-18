args = commandArgs(trailingOnly=T)


occs = read.csv(args[1], sep=args[2])

circlesAbsences <- function(obs, rad, count, lat='decimallatitude', lon='decimallongitude'){
  library(dismo)
  circles = circles(obs[,c(lon, lat)], d=rad, lonlat=T)
  pseudoabs = spsample(circles@polygons, count, type='random', iter=100)
  return(pseudoabs@coords)
}

abs = circlesAbsences(occs, 50000, nrow(occs))

write.table(abs, file=paste0(args[1], "_absences.csv"))