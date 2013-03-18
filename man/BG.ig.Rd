\name{BG.ig}
\alias{BG.ig}
\title{
  Get initial guesses for parameters of a homogeneous- or
  non-homogeneous-background model
}
\description{
  This function calculates initial values of the parameters of a
  fluorescence model in order to initialize the \code{nls} regression.
}
\usage{
BG.ig(adu = NULL, phi.type = "all", f.type = "non-param",
      b.type = "unif", pix.unif = NULL) 
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
}
\details{
  The function first computes an estimate of the mean value of the
  \eqn{\phi_i}{\phi_i} (phi.mean), assuming that the mean value of the
  \eqn{f_j}{f_j} is 1. From this value, we get an initial guess for the
  background fluorescence \code{b}, assumed to be identical for all
  pixels. The initial guesses for all \eqn{f_j}{f_j} (normalized
  fluorescence time course f(t), with baseline equal to 1) and the
  \eqn{\phi_i}{\phi_i} are then deduced from the data and the value of
  \code{phi.mean}.
  If \code{f.type} is set to "mono", a linear model (using \code{lm}) is
  fitted to the logarithm of f(t).
  If \code{b.type} is set to "non-unif", the initial guess for the
  background fluorescence of inhomogeneous pixels (the pixels that are
  not in \code{pix.unif}) is set to the initial guess of \code{b}.
}
\value{
  A named list containing the initial guesses for all model parameters.
}
\author{
  Sebastien Joucla \email{joucla_sebastien@yahoo.fr}
}
\seealso{
  \code{\link{BG.fitAll}}, \code{\link{BG.mkFct}}
}
\examples{
## Compute initial guesses for a non-homogeneous background fluorescence
## model, with a normalized fluorescence transient estimated pointwise.
## The parameters are those used for panel A of Figure 3 of the paper.
set.seed(16051981)
pix.wrong <- c(5,10,15,20)
dB.pc <- 0.30
adu.list <- BG.createData(pix.wrong=pix.wrong, dB.pc=dB.pc, nRep=1)
adu <- adu.list[[1]]
tc <- attr(adu, "tc")
tOn <- attr(adu, "tOn")
nPixels <- nrow(adu)

# Model specifications
phi.type <- "all"
f.type <- "non-param"
b.type <- "non-unif"

## Find initial guesses for all model parameters
pix.unif <- (1:nPixels)[-pix.wrong]
IG <- BG.ig(adu = adu, phi.type = "all", f.type = "non-param",
            b.type = "non-unif", pix.unif = (1:nPixels)[-pix.wrong])
}
