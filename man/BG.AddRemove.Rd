\name{BG.AddRemove}
\alias{BG.AddRemove}
\title{
  Determine the pixels of a ROI having a different background fluorescence 
}
\description{
  This low-level function, called iteratively by \code{BG.track}
  performs one step of the Add/Remove procedure in order to update the
  list of pixels that have a different background fluorescence than the
  others, sharing the same background
  fluorescence. Non-homogeneous-background models are fitted to the
  data, starting from the set of "non-homogeneous" pixels passed as
  argument, which results in an updated set of "non-homogeneous"
  pixels.
}
\usage{
BG.AddRemove(adu = NULL, phi.type = "all", f.type = "non-param",
             tc = NULL, tOn = NULL, pix.unif = NULL,
             pix.nonUnif = NULL, nPixels = NULL, p_thresh = 0.025,
             Plot = FALSE, Print = FALSE)
}
\arguments{
  \item{adu}{
    A 2D matrix containing, in each row "j", the fluorescence transient
    recorded in pixel "j" of the ROI.
  }
  \item{phi.type}{
    A string specifying if the \eqn{\phi_i}{\phi_i} are parameters of the
    model (set \code{phi.type} to its default "all" value) or should be
    calculated as a function of the backgound fluorescence and the
    normalized fluorescence transient f(t) (in this case, set
    \code{phi.type} to "fct_b").
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
  \item{pix.unif}{
    Vector of the indices of the pixels sharing the same (homogeneous)
    background fluorescence, before the application of the function
    (corresponding to the next step of the Add/Remove procedure).
  }
  \item{pix.nonUnif}{
    Vector of the indices of the pixels sharing the same (homogeneous)
    background fluorescence, before the application of the function
    (corresponding to the next step of the Add/Remove procedure).
  }
  \item{nPixels}{
    Number of pixels of the Region Of Interest.
  }
  \item{p_thresh}{
    A numerical value (strictly below 1) corresponding to the threshold
    p-value (default value 0.025) under which pixels tested for
    background inhomogeneity are considered to have a background
    fluorescence statistically different from the common one.
  }
  \item{Plot}{
    A string specifying in which type of device the results of the
    Add/Remove procedure should be plotted. This can be "x11", "pdf",
    "jpeg" or "tiff". Using any other value, no device will be opened.
    Note that if \code{Plot } is set to "pdf", "jpeg" or "tiff", the
    device will be closed at the end of the procedure.
  }
  \item{Print}{
    A logical value indicating if the pixels tested for background
    inhomogeneity should be printed in the terminal.
  }
}
\details{
  All details of the Add/Remove procedure can be found in the article.
}
\value{
  A vector of pixels (between 1 and \code{nPixels}) having a
  statistically different background fluorescence.
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
  \code{\link{BG.track}}
}
\examples{
## This example reproduces the first step of the analysis leading to
## Figure 5 of the paper (Joucla et al., Cell Calcium, in revision)

## Create a single data set
## (fluorescence parameters are taken as the default arguments values)
set.seed(16051981)
pix.wrong <- c(5,10,15,20)
dB.pc <- 0.30
adu.list <- BG.createData(pix.wrong=pix.wrong, dB.pc=dB.pc, nRep=1)
adu <- adu.list[[1]]
tc <- attr(adu, "tc")
tOn <- attr(adu, "tOn")

## Fit the whole model on all data (uniform background)
directFit <- BG.fitAll(adu, phi.type="all", b.type="unif",
                       pix.unif=NULL, f.type="mono", tc, tOn)

## Step1: Start with all pixels homogeneous except 1
## compute a p-value for each of them
nPixels <- nrow(adu)
pix.unif <- 1:nPixels
pix.unif.bk <- pix.unif
pix.nonUnif <- c()
layout(matrix(1:4, nrow=1))
pix2exclude <- BG.AddRemove(adu, phi.type="all", f.type="mono",
                            tc, tOn, pix.unif, pix.nonUnif, nPixels,
                            p_thresh=0.025, Plot="x11", Print=TRUE)
}
