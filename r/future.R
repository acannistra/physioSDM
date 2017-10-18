## Physiology in Species Distribution Models
## Tony Cannistra & Lauren Buckley, 2017
## ~ Future Projection and Evaluation helper functions. 
suppressPackageStartupMessages({
  library(raster)
  library(plyr)
})

getClumps <- function(predictionRaster, threshold=0.5){
  predictionRaster[predictionRaster >= threshold] = 1
  predictionRaster[predictionRaster < threshold] = NA
  return(clump(predictionRaster, directions=4))
}

largestClumpId <- function(clumps, threshold=20){
  #get number of cells in eeach clump
  freqs <- as.data.frame(freq(clumps, useNA = 'no'))
  #return clump ID of largest clump
  return(freqs[which(freqs$count == max(freqs$count)),]$value)
}
  
plotClumps <- function(rc){
  clump_id <- getValues(rc)    
  xy <- xyFromCell(rc,1:ncell(rc))
  df <- data.frame(xy, clump_id, is_clump = rc[] %in% freq(rc, useNA = 'no')[,1])
  df[df$is_clump == T, ]
  dfm <- ddply(df[df$is_clump == T, ], .(clump_id), summarise, xm = mean(x), ym = mean(y))
  plot(rc)
  text(dfm[, 2:3], labels = dfm$clump_id)
}

comparePredictions <- function(predictions, threshold=0.5, titles=NA, occs=NA, verbose=FALSE, 
                               fullresult=FALSE, arrows=TRUE){
  npreds = length(predictions)
  par(mar=c(5, 4.6, 4, 2.6))
  #plot.new()
  orig_x = NA
  orig_y = NA
  if(fullresult){
    par(mfrow=c(2, npreds))
  } else{
    par(mfrow=c(1, npreds))
  }
  if(fullresult){
    for (p in seq_along(predictions)){
      if(is.na(threshold)){
        plot(predictions[[p]], 
             legend.width=1.9, 
             cex.axis=2,
             cex.lab=1.6,
             xlab='Longitude (deg)', ylab="Latitude(deg)")
      } else {
        plot(predictions[[p]] >= threshold,
             legend.width=1.9,
             cex.axis=2,
             cex.lab = 1.6,
             xlab='Longitude (deg)', ylab="Latitude(deg)")
      }
      if(!is.na(titles)){
        title(titles[p], cex.main = 1.9)
      }
      if(!is.na(occs)){
        points(x=occs$x, y=occs$y, pch='.')
      }
    }
  }
  p_idx = 0
  for(p in predictions){
    p_idx = p_idx + 1
    pclumps = getClumps(p)
    largestclump = tryCatch((pclumps == largestClumpId(pclumps)), 
                            error=function(e){
                              flog.error(sprintf("failure in circles: %s", str(e$message)))
                              return(NULL)
                            })
    if(is.null(largestclump)){
      return()
    }
    plot(largestclump, 
         cex = 1.5, 
         cex.axis = 2,
         cex.lab = 1.6,
         legend=F,
         xlab='Longitude (deg)', ylab="Latitude(deg)")
    pcells = as.data.frame(rasterToPoints(largestclump))
    pcells = pcells[pcells$layer == 1, ]
    center_x = mean(pcells$x)
    center_y = mean(pcells$y)
    if(!fullresult){
      if(!is.na(titles)){
        title(titles[p_idx], cex.main = 1.9)
      }
    }
    if(verbose){
      title(sprintf("# cells: %f\n(%f, %f)", nrow(pcells),center_x, center_y))
    }
    if(is.na(orig_x) && is.na(orig_y)){
      orig_x = center_x
      orig_y = center_y
    } else{
      if(arrows){
        points(x=orig_x, y=orig_y, col='red')
        arrows(x0=orig_x, y0=orig_y, x1=center_x, y1 = center_y, length=0.1)
      }
    }
    points(x=mean(pcells$x), y=mean(pcells$y))
    if(!is.na(occs)){
      points(x=occs$x, y=occs$y, pch='.')
    }
    
  }
}
