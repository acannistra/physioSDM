library(maps)
library(mapdata)
library(ggmap)
library(dismo)
library(dplyr)


circlesAbsences <- function(obs, rad, count, lat='decimallatitude', lon='decimallongitude'){
  library(dismo)
  circles = circles(obs[,c(lon, lat)], d=rad, lonlat=T)
  pseudoabs = spsample(circles@polygons, count, type='random', iter=100)
  return(pseudoabs@coords)
}

appendDates <-function(pts, obs, lat='decimallatitude', lon='decimallongitude'){
  obs_pts = SpatialPoints(cbind(obs[,lon], obs[,lat]))
  tree = createTree(obs_pts)
  og_inds = knnLookup(tree, newdat = pts, k=1)
  pts = data.frame(pts) %>% 
        mutate(month=obs[og_inds,'month'], day=obs[og_inds,'day'], year=obs[og_inds,'year'])
  return(pts)
}

occs_all = read.csv('../data/occs/sceloporus_occidentalis_westus.csv', sep='\t')

filtered = occs_all %>% filter(coordinateuncertaintyinmeters < 5000) %>%
                       filter(!((year < 2000) & (coordinateuncertaintyinmeters < 3))) %>%
                       filter(!((year >= 2000) & (coordinateuncertaintyinmeters < 200)))

early = filtered %>% filter((year <=1939) & (year >=1900))
late = filtered %>% filter((year >= 1970) & (year <= 2009))

early_abs = circlesAbsences(early, 50000, nrow(early))




#ggmap(get_stamenmap(extent, maptype='terrain-lines', zoom=6)) +
#      geom_point(aes(x=x, y=y), data=as.data.frame(pseudoabs@coords), color='darkred') + 
#      geom_point(aes(x=decimallongitude, y=decimallatitude), data=highres, color='blue')
