\name{BG.plotAllPixels}
\alias{BG.plotAllPixels}
\title{
  Plot the raw and predicted fluorescence transients of each individual
  pixel using a 2D representation.
}
\description{
  This function plots the individual fluorescence transients recorded in
  each pixel, using a 2D representation. The data ("raw" data in black,
  "predicted" data in a user-defined color) are plotted in a single
  screen (no \code{layout} or \code{split.screen} is used is used)
}
\usage{
BG.plotAllPixels(BgNonUnifFitObject, ROI, what = "predicted", color = "red")
}
\arguments{
  \item{BgNonUnifFitObject}{
    A \code{nls} object returned by either the \code{BG.fitAll} or
    \code{BG.track} function.
  }
  \item{ROI}{
    A list of two elements
    - \code{i} a vector of integers corresponding to the "x" coordinates
    of the pixels used for the fit (takne from the original image data).
    - \code{j} a vector of integers corresponding to the "y" coordinates
    of the pixels used for the fit (taken from the original image data).
  }
  \item{what}{
    A string specifying if "raw" data should be plotted alone, or if
    "predicted" data should be superimposed.
  }
  \item{color}{
    A string specifying which color to use (useful only if \code{"what"})
    is set to "predicted".
  }
}
\value{
  No value is returned.
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
\examples{
## Reproduce panel B of Figure 9 of the paper
rm(list=ls())

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

BG.plotAllPixels(directFit_h, POMC$ROI, "raw", "red")
}