\name{BG.createData}
\alias{BG.createData}
\title{
Create data sets according to the fluorescence model described in Joucla
et al. (Cell Calcium, in revision).
}
\description{
This function creates several sets of fluorescence data according to the
fluorescence model described in Joucla et al. (Cell Calcium, in
revision). This model embeds a normalized fluorescence time course f(t)
modeled with a monoexponential return to baseline. The data are drawn
from Poisson distribution with parameters equal to the true fluorescence
value (for each pixel in the ROI and each time sample) and are then
multiplied by the gain of the camera.
}
\usage{
BG.createData(nPixels = 20, phi = seq(0, 1, length.out = 20 + 1)[-1],
Factor = 500, F0 = 1, dF = 0.5, tau = 1, tOn = 1, tc = seq(0, 5, 0.1), b
= 250, dB.pc = 0, pix.wrong = NULL, G = 0.146, nRep = 1)
}
\arguments{
  \item{nPixels}{
Number of pixels of the Region Of Interest (integer). The pixels are
then indexed from 1 to \code{nPixels}.
}
  \item{phi}{
Vector of numerical values corresponding to the amplitude of the
pixel-specific scaling factors. These values are normalized between 0
and 1.
}
  \item{Factor}{
This numerical value is the multiplying factor to apply to vector
\code{phi}, so as to obtain absolute (and not normalized) values.
}
  \item{F0}{
Fluorescence value at the baseline. Generally set to 1, so that the
fluorescence time course f(t) is directly read as relative fluorescence
values.
}
  \item{dF}{
Amplitude of the fluorescence rise at tc=tOn, compared to a baseline
fluorescence equal to 1.
}
  \item{tau}{
Time constant of the fluorescence return to baseline (in second).
}
  \item{tOn}{
Time of the fluorescence rise (in second).
}
  \item{tc}{
Vector of time at which the fluorescence values are observed (in
seconds).
}
  \item{b}{
Amplitude of the background fluorescence. This value is set to all
pixels of the ROI, except those in vector \code{pix.wrong} (in the
latter, the background fluorescence is given by b*(1+dB.pc)).
}
  \item{dB.pc}{
Background fluorescence variation, in percentage of the background
fluorescence common to most pixels of the ROI.
}
  \item{pix.wrong}{
Vector of indices of the pixels of the ROI that have a different
background fluorescence.
}
  \item{G}{
Gain of the acquisition system (a CCD camera for instance).
}
  \item{nRep}{
Number of data sets to create.}
}
\details{
For each data set, the ideal (noise-free) fluorescence time course
\eqn{F_i(t)}{F_i(t)} in pixel i as:
\deqn{%
F_i(t) = \phi_i \cdot f(t) + b
}{%
F_i(t) = phi_i * f(t) + b
},
where \eqn{\phi_i}{phi_i} is the \eqn{i^{th}}{i^{th}} value in vector
\eqn{\phi}{phi}.
The above expression stands for pixels with homogeneous background
\eqn{b}{b}. For "non-homogeneous" pixels (those in vector
\code{pix.wrong}), the background fluorescence is \eqn{b \cdot
  (1+\code{dB.pc})}{b*(1+dB.pc)}.
The "noisy" data are then generated according to a Poisson distribution
and multiplied by the gain \eqn{G}{G} of the acquisition system.
}
\value{
A list of \code{nRep} elements, each of them being a matrix with
\code{nPixels} rows and length(tc) columns. Each element has the
following attributes:
  \item{\code{param.true}}{A named vector giving the values of the
    parameters used to create the data}
  \item{\code{f}}{A vector corresponding to the normalized fluorescence
    time course f(t)}
  \item{\code{tc}}{A copy of argument \code{tc}}
  \item{\code{tOn}}{A copy of argument \code{tOn}}
  \item{\code{nPixels}}{A copy of argument \code{nPixels}}
}
\author{
Sebastien Joucla \email{joucla_sebastien@yahoo.fr}
}
\examples{
## This example reproduces Figure 1 of the paper (Joucla et al., Cell
## Calcium, in revision)

## Create a single data set
## (fluorescence parameters are taken as the default arguments values)
set.seed(16051981)
adu.list <- BG.createData(nRep=1)

## Get data from adu.list
adu <- adu.list[[1]]
f <- attr(adu, "f")
tc <- attr(adu, "tc")
tOn <- attr(adu, "tOn")
P <- attr(adu, "param.true")
b <- P["b"]
Struct <- strsplit(names(P), "phi_")
phi.idx <- which(is.finite(sapply(1:length(Struct),
                                  function(i)
                                    as.numeric(Struct[[i]][2]))))
phi <- P[phi.idx]

## Reorganize the data according to a 4*5 ROI
phi.mat <- matrix(phi, nrow=4, byrow=TRUE)
phi.mat <- phi.mat[seq(nrow(phi.mat),1,-1),]

b.mat <- matrix(b, nrow=4, ncol=5, byrow=TRUE)

x.image <- 1:5
y.image <- 1:4
x.rep <- matrix(data = rep(x.image, length(y.image)),
                nrow=length(y.image), byrow=TRUE)
y.rep <- matrix(data = rep(y.image, each=length(x.image)),
                nrow=length(y.image), byrow=TRUE)

adu.array <- array(NA, dim=c(4,5,51))
for(j in 1:4) {
  for(i in 1:5) {
    adu.array[j,i,] <- adu[i+5*(j-1),]
  }
}

## Open the window to create Figure 1
X11(width=7, height=3.4)
mat <- rbind(c(1,2,5),c(3,3,5),c(4,4,5))
layout(mat, widths = c(1.5,1.5,4), heights=c(1.5,0.42,1.5))
par(oma=c(0.5,0,1,0.5))

## A : Background image
par(mar=c(0,1,1.5,1))
image(x.image, y.image, t(b.mat)/max(phi.mat), c(0,1),
      col=gray.colors(21,0,1), axes=FALSE, ann=FALSE)
for(idx.x in 1:dim(x.rep)[1]) {
  for(idx.y in 1:dim(y.rep)[2]) {
    text(x.rep[idx.x, idx.y]+0.6,
         y.rep[max(y.rep)-idx.x+1, idx.y]+0.3,
         idx.y + (idx.x-1)*dim(x.rep)[2], pos=2, col=rgb(0.4,0.4,0.4), cex=0.75)
  }
}
mtext("   Background image", side=3, line=0.5, cex=0.75)
mtext("A", side=3, at=0.6, line=0.5, font=2, cex=0.9)
abline(v=x.image[1:(length(x.image)-1)]+0.5, lty=2, lwd=0.3)
abline(h=y.image[1:(length(y.image)-1)]+0.5, lty=2, lwd=0.3)
box()

## B : Image of the phi_i(s)
par(mar=c(0,1,1.5,1))
image(x.image, y.image, t(phi.mat)/max(phi.mat), c(0,1),
      col=gray.colors(21,0,1), axes=FALSE, ann=FALSE)
for(idx.x in 1:dim(x.rep)[1]) {
  for(idx.y in 1:dim(y.rep)[2]) {
    text(x.rep[idx.x, idx.y]+0.6,
         y.rep[max(y.rep)-idx.x+1, idx.y]+0.3,
         idx.y + (idx.x-1)*dim(x.rep)[2], pos=2,
         col=rgb(0.4,0.4,0.4), cex=0.75)
  }
}
mtext(expression(paste("Image of the ", phi[i])),
      side=3, line=0.5, cex=0.75)
mtext("B", side=3, at=0.6, line=0.5, font=2, cex=1)
abline(v=x.image[1:(length(x.image)-1)]+0.5, lty=2, lwd=0.3)
abline(h=y.image[1:(length(y.image)-1)]+0.5, lty=2, lwd=0.3)
box()

## Colorbar
par(mar=c(1.25,1,1,1))
x_color <- 0:length(b.mat)
image(x_color, 1, as.matrix(sort(x_color, decreasing=FALSE))/max(x_color),
      c(0,1), col=gray.colors(length(x_color)+1,0,1), axes=FALSE, ann=FALSE)
axis(1, at=c(min(x_color), median(x_color), max(x_color)),
     labels=c(min(x_color), median(x_color), max(x_color))*max(phi.mat)/max(x_color),
     tick=FALSE, mgp=c(3,0,0), cex.axis=0.8)
box()

## C : Plot of the normalized time course
par(mar=c(1.5,2.25,2.25,0.5))
plot(tc, f, type="l", axes=FALSE, ann=FALSE, xaxs="i", xlim=c(0,5.25))
axis(1, at=0:5, mgp=c(3,0.5,0), cex.axis=0.8)
axis(2, at=seq(min(f), max(f), by=0.1), mgp=c(3,0.75,0),
     labels=seq(min(f), max(f), by=0.1), las=1, cex.axis=0.8)
axis(1, at=5.25, col.ticks="white", col="white", labels="(s)",
     mgp=c(3,0.5,0), cex.axis=0.8)
mtext("Normalized fluorescence time course f(t)",
      side=3, line=0.5, cex=0.75)
mtext("C", side=3, at=-0.2, line=0.5, font=2, cex=0.9)

## D : Plot of the fluorescence time course in all pixels
xLim <- c(min(x.image)-0.5, max(x.image)+0.5)
yLim <- c(min(y.image)-0.5, max(y.image)+0.5)

par(mar=c(0.25,3.5,1.5,2.25))
plot(1, 1, type="n", axes=FALSE, ann=FALSE, xlim = xLim,
     ylim = yLim, xaxs="i", yaxs="i")

factor <- 1250

for(idx.x in 1:dim(x.rep)[1]) {
  xx <- NULL
  yy <- NULL
  for(idx.y in 1:dim(y.rep)[2]) {
    x_ref <- x.rep[idx.x, idx.y]
    ## print(x_ref)
    y_ref <- y.rep[max(y.rep)-idx.x+1, idx.y]
    ## print(y_ref)
    xx <- c(xx, NA, x_ref +
            (1:dim(adu.array)[3])/dim(adu.array)[3] - 0.5)
    yy <- c(yy, NA, y_ref +
            adu.array[y.rep[idx.x, idx.y], x.rep[idx.x, idx.y],]/factor-0.5)
    text(x.rep[idx.x, idx.y]+0.55,
         y.rep[max(y.rep)-idx.x+1, idx.y]+0.4,
         idx.y + (idx.x-1)*dim(x.rep)[2], pos=2, col=rgb(0.4,0.4,0.4))
    ## print(idx.y + (idx.x-1)*dim(x.rep)[2])
  }
  lines(xx, yy, col="black")
}
axis(4, at=seq(0.5,1.3, length.out=5), labels=seq(0, factor-1, by=250),
     mgp=c(3,0.75,0), las=1, cex.axis=0.8)
axis(4, at=2+seq(0.5,1.3, length.out=5), labels=seq(0, factor-1, by=250),
     mgp=c(3,0.75,0), las=1, cex.axis=0.8)
axis(2, at=1+seq(0.5,1.3, length.out=5), labels=seq(0, factor-1, by=250),
     mgp=c(3,0.75,0), las=1, cex.axis=0.8)
axis(2, at=3+seq(0.5,1.3, length.out=5), labels=seq(0, factor-1, by=250),
     mgp=c(3,0.75,0), las=1, cex.axis=0.8)

mtext(expression(paste("Raw fluorescence time courses ", F[i], "(t)")),
      side=3, line=0.5, cex=0.75)
mtext("D", side=3, at=0.3, line=0.5, font=2, cex=0.9)

## Pixels limits
abline(v = seq(xLim[1]+1, xLim[2]-1, 1), col="gray")
abline(h = seq(yLim[1]+1, yLim[2]-1, 1), col="gray")

## Figure borders
box()
}
