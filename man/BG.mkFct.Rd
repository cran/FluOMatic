\name{BG.mkFct}
\alias{BG.mkFct}
\title{
Create the function passed to nls in order to fit the fluorescence model.
}
\description{
  This low-level function, called by \code{BG.fitAll} and
  \code{BG.track}, creates a function that generates fluorescence data
  according to the current fluorescence model (described in
  \code{\link{BG.createData}}). It is then passed to \code{nls} to
  perform the nonlinear regression. 
}
\usage{
BG.mkFct(adu = NULL, phi.type = "all", b.type = "unif", pix.unif = NULL,
         f.type = NULL, SQRT = TRUE)
}
\arguments{
  \item{adu}{
    A 2D matrix containing, in each row "j", the fluorescence transient
    recorded in pixel "j" of the ROI.
  }
  \item{phi.type}{
    A string specifying if the \eqn{phi_i}{phi_i} are parameters of the
    model (set \code{phi.type} to its default "all" value) or should be
    calculated as a function of the backgound fluorescence and the
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
  \item{SQRT}{
    A logical value specifying wether the square-root transformation
    should be used (set \code{SQRT} to TRUE, the default value) or not
    (set \code{SQRT} to FALSE).
  }
}
\value{
  A function with arguments corresponding to the model parameters. The
  latter depend on the type of background fluorescence ("unif" or
  "non-unif"), normalized fluorescence transient ("mono" or
  "non-param") and pixel-specific amplitude factors ("all" or "fct_b").
}
\author{
  Sebastien Joucla \email{joucla_sebastien@yahoo.fr}
}
\seealso{
  \code{\link{BG.ig}}
}
\examples{
## Create a function to pass to nls in order to perform a non-linear
## regression on fluorescence data
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

SQRT <- TRUE
fitFct <- BG.mkFct(adu, phi.type, b.type, pix.unif, f.type, SQRT)
}
