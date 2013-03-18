\name{FluOMatic-package}
\alias{FluOMatic-package}
\alias{FluOMatic}
\docType{package}
\title{Estimation of background-subtracted fluorescence data}
\description{Simulate and analyse fluorescence data sets obtained from
  calcium imaging experiments in living cells (such as neurons). The
  package provide tools to estimate background-subtracted fluorescence
  transients and build meaningful confidence intervals on these
  estimations. This covers both the case of homogeneous background
  fluorescence and that of Regions Of Interest (ROI) displaying
  background inhomogeneities. Background-subtracted fluorescence
  transients can be fitted with monoexponential models or
  estimated pointwise.}
\details{
\tabular{ll}{
Package: \tab FluOMatic\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-03-18\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Sebastien Joucla, Christophe Pouzat

Maintainer: Sebastien Joucla <joucla_sebastien@yahoo.fr>
}
\references{
Joucla S, Franconville R, Pippow A, Koppenburg P, Pouzat C (2013)
Estimating background-subtracted fluorescence transients in calcium
imaging experiments: a quantitative approach (Cell Calcium, in
revision)}
\keyword{ package }
\seealso{
  see package \code{CalciOMatic}.
}
\examples{
## Reproduce the analysis leading to Figure 9 of the paper

## Load the POMC data set
data(POMC)

## Build the (2D) adu matrix from the pixels in the selected ROI
nPixels <- length(POMC$ROI$i)
adu <- t(sapply(1:length(POMC$ROI$i),
                function(k)
                POMC$adu.img[POMC$ROI$i[k], POMC$ROI$j[k], ]))

## Fit a homogeneous-background model, using a non-parametric
## description of f(t) (pointwise estimation)
phi.type <- "all"
b.type <- "unif"
pix.unif <- 1:length(POMC$ROI$i)
f.type <- "non-param"
directFit_h <- BG.fitAll(adu, phi.type, b.type, pix.unif, f.type,
                         POMC$time, POMC$tOn)
attr(directFit_h, "RawData") <- adu

\dontrun{
## Fit a non-homogeneous-background model, using a non-parametric
## description of f(t) (pointwise estimation)
## Define the threshold p-value for the exclusion of non-homogeneous pixels
p_thresh <- 0.025
## Using "x11" as last argument, the list of pixels having a background fluorescence
## different from that of most of the pixels of the ROI will be displayed at each step
## of the "Add/Remove" algorithm
directFit_nh <- BG.track(adu, f.type, phi.type, pix.unif, POMC$time,
                         POMC$tOn, "x11", TRUE)
attr(directFit_nh, "RawData") <- adu
## List of pixels having a statistically non-different background fluorescence
pix.unif <- attr(directFit_nh, "pix.unif")
}
}
