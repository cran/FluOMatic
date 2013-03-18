BG.plotAllPixels <-
function(BgNonUnifFitObject, ROI, what="predicted", color="red") {
  ## Function plot.BgNonUnifFitObject
  
  ## Get the raw and predicted data, as well as the residuals
  adu.raw <- attr(BgNonUnifFitObject, "RawData")
  if(is.null(adu.raw)) adu.raw <- attr(BgNonUnifFitObject, "adu")
  adu.fit <- matrix(predict(BgNonUnifFitObject), nrow(adu.raw))
  
  ## Open the figure to plot either the raw and predicted data, or the residuals
  ## The pixel time courses will be plotted in a single layout
  
  ## Prepare the pixels coordinates
  x.image <- min(ROI$i):max(ROI$i)
  y.image <- min(ROI$j):max(ROI$j)
  xLim <- c(min(x.image)-0.5, max(x.image)+0.5)
  yLim <- c(min(y.image)-0.5, max(y.image)+0.5)
  plot(1, 1, type="n", xlim = xLim, ylim = yLim, xaxs="i", yaxs="i", axes=FALSE, ann=FALSE)
  
  ## Prepare the data
  nPixels <- length(ROI$i)
  xx <- c()
  yy.raw <- c()
  yy.fit <- c()
  for(k in 1:nPixels) {
    ## The x data
    xx <- c(xx, NA, ROI$i[k]+(1:ncol(adu.raw))/ncol(adu.raw)-0.5)
    ## The raw y data
    y.toadd <- sqrt(adu.raw[k,])
    mini <- min(y.toadd)
    y.toadd <- y.toadd-mini
    maxi <- max(y.toadd)
    y.toadd <- y.toadd/maxi
    y.toadd <- 0.04+0.92*y.toadd
    yy.raw <- c(yy.raw, NA, ROI$j[k]+y.toadd-0.5)
    ## The fitted y data
    y.toadd <- adu.fit[k,]
    y.toadd <- y.toadd-mini
    y.toadd <- y.toadd/maxi
    y.toadd <- 0.04+0.92*y.toadd
    yy.fit <- c(yy.fit, NA, ROI$j[k]+y.toadd-0.5)
  }
  
  ## Plot the data
  if(what == "predicted" | what == "raw") {
    lines(xx, yy.raw, col="black")
  }
  if(what == "predicted") {
    lines(xx, yy.fit, col=color)
  }
  ## Add a box around each pixel
  for(k in 1:nPixels) {
    rect(ROI$i[k]-0.5, ROI$j[k]-0.5,ROI$i[k]+0.5, ROI$j[k]+0.5, border="gray")
  }
  
  ## Pixels limits
  ## abline(v = seq(xLim[1]+1, xLim[2]-1, 1), col="gray")
  ## abline(h = seq(yLim[1]+1, yLim[2]-1, 1), col="gray")
    
  ## Figure borders
#    axis(side=1, at=seq(xLim[1]+0.5, xLim[2]-0.5, 1), las=1, tcl=0,
#         mgp=c(3, 0.25, 0), labels=unique(as.vector(x.rep)))
#   axis(side=2, at=seq(yLim[1]+0.5, yLim[2]-0.5, 1), las=1, tcl=0,
#         mgp=c(3, 0.25, 0), labels=unique(as.vector(y.rep)))
  box()
}

