\name{BG.fitAll}
\alias{BG.fitAll}
\title{
Fit a homogeneous- or non-homogeneous-background model on a single-trial
data set
}
\description{
This function fits a homogeneous- or non-homogeneous-background model on
a fluorescence transient in order to estimate the background
fluorescence, the pixel-sopecific amplitude factors and the normalized
fluorescence time course f(t). In the case of a non-homogeneous
background fluorescence, the indices of the inhomogeneous pixels have to
be specified.
}
\usage{
BG.fitAll(adu = NULL, phi.type = "all", b.type = "unif", pix.unif =
NULL, f.type = NULL, tc = NULL, tOn = NULL)
}
\arguments{
  \item{adu}{
    A 2D matrix containing, in each row "j", the fluorescence transient
    recorded in pixel "j" of the ROI.
  }
  \item{phi.type}{
    A string specifying if the \eqn{\phi_i}{\phi_i} are parameters of
    the model (set \code{phi.type} to its default "all" value) or should
    be calculated as a function of the backgound fluorescence and the
    normalized fluorescence transient f(t) (in this case, set
    \code{phi.type} to "fct_b").
  }
  \item{b.type}{
    A string specifying the type of background fluorescence to consider:
    either uniform (set \code{b.type} to "unif") or non-uniform (set
    \code{b.type} to "non-unif"). In the first case, argument
    \code{pix.unif} is ignored. In the second case, the indices of the
    pixels sharing the same background fluorescence should be specified in
    argument \code{pix.unif}.
  }
  \item{pix.unif}{
    Vector of the indices of the pixels sharing the same (homogeneous)
    background fluorescence (has to be specified only if \code{b.type}
    is "non-unif").
  }
  \item{f.type}{
    A string specifying if the normalized fluorescence transient f(t) is
    modeled by a monoexponential model (set \code{f.type} to "mono") or
    is estimated pointwise (set \code{f.type} to "non-param").
  }
  \item{tc}{
    A vector of times at which the fluorescence data are observed
    (useful only in the case where \code{f.type} is set to "mono").
  }
  \item{tOn}{
    A numerical value corresponding to the time at which the
    fluorescence instantanously rises from the baseline to its maximal
    value (useful only in the case where \code{f.type} is set to
    "mono").
  }
}
\details{
  The fluorescence model fitted to the data corresponds to that
  described in \code{\link{BG.createData}}.
}
\value{
  A \code{nls} object containing all the informations of the non-linear
  regression.
}
\references{
  Joucla S, Franconville R, Pippow A, Koppenburg P, Pouzat C (2013)
  Estimating background-subtracted fluorescence transients in calcium
  imaging experiments: a quantitative approach (Cell Calcium, in
  revision)
}
\author{
  Sebastien Joucla \email{joucla_sebastien@yahoo.fr}
}
\seealso{
  \code{\link{BG.createData}}, \code{\link{BG.track}}
}
\examples{
## This example reproduces panel A of Figure 2 of the paper
## (Joucla et al., Cell Calcium, in revision)

## Create a single data set
## (fluorescence parameters are taken as the default arguments values)
set.seed(16051981)
adu.list <- BG.createData(nRep=1)
adu.backup <- adu.list[[1]]
tc <- attr(adu.backup, "tc")
tOn <- attr(adu.backup, "tOn")
nPixels <- nrow(adu.backup)

# Model specifications
phi.type <- "all"
f.type <- "mono"

## Fit the whole model on all data (uniform background)
b.type <- "unif"
pix.unif <- 1:nPixels
directFit <- BG.fitAll(adu.backup, phi.type, b.type, pix.unif, f.type, tc, tOn)
attr(directFit, "adu") <- adu.backup

## ===================
## Retrieve fitted adu
## ===================

## Get all fitted parameters
fitted.param <- summary(directFit)$parameters[,1]

## Choice of the fluorescence transient f(t), normalize it
if(f.type == "mono") {
  F0 <- 1
  fitted.dF <- fitted.param["dF"]
  fitted.tau <- fitted.param["tau"]
  fitted.f <- F0 + fitted.dF * ifelse(tc >= tOn, exp(-(tc-tOn)/fitted.tau), 0)
} else if(f.type == "non-param") {
  Struct <- strsplit(names(fitted.param), "f_")
  f.idx <- sapply(1:length(Struct), function(i) as.numeric(Struct[[i]][2]))
  fitted.f <- c(1,unlist(sapply(2:(length(which(is.finite(f.idx)))+1),
                                function(i)
                                  fitted.param[which(f.idx == i)])))
}

## Background fluorescence
fitted.b <- fitted.param["b"]
fitted.b <- rep(fitted.b, nPixels)

## Background fluorescence on the non-uniform pixels
Struct <- strsplit(names(fitted.param), "b_")
b.idx <- sapply(1:length(Struct), function(i) as.numeric(Struct[[i]][2]))
b.idx2 <- b.idx[which(is.finite(b.idx))]
b.idx <- which(is.finite(b.idx))
fitted.b[b.idx2] <- fitted.param[b.idx]

## Pixel-specific scaling factors phi_i
Struct <- strsplit(names(fitted.param), "phi_")
phi.idx <- sapply(1:length(Struct), function(i) as.numeric(Struct[[i]][2]))
fitted.phi <- sapply(1:length(which(is.finite(phi.idx))),
                     function(i)
                       fitted.param[which(phi.idx == i)])

## Construction of the whole ideal adu matrix
fitted.adu <- fitted.phi %o% fitted.f
fitted.B <- matrix(fitted.b, nrow(adu.backup), ncol(adu.backup), byrow=!TRUE)

## Add the background to get the complete unnoised data
fitted.adu <- fitted.adu + fitted.B

## ======================================
## Plot the fitted adu above the raw data
## ======================================

## Prepare the data to plot
x.image <- 1:5
y.image <- 1:4
x.rep <- matrix(data = rep(x.image, length(y.image)),
                nrow=length(y.image), byrow=TRUE)
y.rep <- matrix(data = rep(y.image, each=length(x.image)),
                nrow=length(y.image), byrow=TRUE)

adu.array <- array(NA, dim=c(4,5,length(tc)))
for(j in 1:4) {
  for(i in 1:5) {
    adu.array[j,i,] <- adu.backup[i+5*(j-1),]
  }
}

fitted.adu.array <- array(NA, dim=c(4,5,length(tc)))
for(j in 1:4) {
  for(i in 1:5) {
    fitted.adu.array[j,i,] <- fitted.adu[i+5*(j-1),]
  }
}

## Open the window
## X11(width=3.5, height=3.5)
## par(oma=rep(0.5,4))
## xLim <- c(min(x.image)-0.5, max(x.image)+0.5)
## yLim <- c(min(y.image)-0.5, max(y.image)+0.5)
xLim <- c(1.5,4.5)
yLim <- c(0.5,3.5)

par(mar=c(0.5,0,1,0.5))
plot(1, 1, type="n", axes=FALSE, ann=FALSE, xlim = xLim,
     ylim = yLim, xaxs="i", yaxs="i")

offset <- 350
factor <- 650
step <- 250

for(idx.x in 2:4) {
  xx <- NULL
  yy <- NULL
  fitted.yy <- NULL
  for(idx.y in 2:4) {
    x_ref <- x.rep[idx.x, idx.y]
    y_ref <- y.rep[max(y.rep)-idx.x+1, idx.y]
    xx <- c(xx, NA, x_ref +
      (1:dim(adu.array)[3])/dim(adu.array)[3] - 0.5)
    yy <- c(yy, NA, y_ref +
      (adu.array[y.rep[idx.x, idx.y], x.rep[idx.x, idx.y],]-offset)/factor-0.5)
    fitted.yy <- c(fitted.yy, NA, y_ref +
      (fitted.adu.array[y.rep[idx.x, idx.y], x.rep[idx.x, idx.y],]-offset)/factor-0.5)
    text(x.rep[idx.x, idx.y]+0.57,
         y.rep[max(y.rep)-idx.x+1, idx.y]+0.4,
         idx.y + (idx.x-1)*dim(x.rep)[2], pos=2, col=rgb(0.4,0.4,0.4), cex=0.7)
  }
  lines(xx, yy, col="black")
  lines(xx, fitted.yy, col="red")
}

axis(2, at=0+seq(0.5,1.45, by=step/factor), labels=seq(offset, offset+factor-1, by=step),
     mgp=c(3,0.7,0), las=1, cex.axis=0.7)
axis(2, at=1+seq(0.5,1.45, by=step/factor), labels=seq(offset, offset+factor-1, by=step),
     mgp=c(3,0.7,0), las=1, cex.axis=0.7)
axis(2, at=2+seq(0.5,1.45, by=step/factor), labels=seq(offset, offset+factor-1, by=step),
     mgp=c(3,0.7,0), las=1, cex.axis=0.7)

mtext("Homogeneous background data", side=3, line=0.25, cex=0.9)
mtext("A", side=3, at=0.94, line=0.25, font=2, cex=0.9)
mtext("Raw (black) and fitted (red) time courses", side=2, line=2.25, cex=0.8)

## Pixels limits
abline(v = seq(xLim[1]+1, xLim[2]-1, 1), col="gray")
abline(h = seq(yLim[1]+1, yLim[2]-1, 1), col="gray")

## Figure borders
box()
}
