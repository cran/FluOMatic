\name{BG.track}
\alias{BG.track}
\title{
  Fit the data with a non-homogeneous-background fluorescence model
}
\description{
  This function fits data set with a non-uniform background model. It is
  based on a recursive approach aiming at finding the pixels that have a
  background fluorescence different than the others.
}
\usage{
BG.track(adu = NULL, f.type = "mono", phi.type = "all", pix.unif = NULL,
         tc = NULL, tOn = NULL, p_thresh = 0.025, Plot = "x11",
         Print = FALSE)
}
\arguments{
  \item{adu}{
    A 2D matrix containing, in each row "j", the fluorescence transient
    recorded in piel "j" of the ROI.
  }
  \item{f.type}{
    A string specifying if the normalized fluorescence transient f(t) is
    modeled by a monoexponential model (set \code{f.type} to "mono") or
    is estimated pointwise (set \code{f.type} to "non-param").
  }
  \item{phi.type}{
    A string specifying if the \eqn{phi_i}{phi_i} are parameters of the
    model (set \code{phi.type} to its default "all" value) or should be
    calculated as a function of the backgound fluorescence and the
    normalized fluorescence transient f(t) (in this case, set
    \code{phi.type} to "fct_b").
  }
  \item{pix.unif}{
    Vector of the indices of the pixels sharing the same (homogeneous)
    background fluorescence. This is generally set to the default NULL
    value, since all pixels are generally assumed to have the same
    background fluorescence at the beginning of the Add/Remove
    procedure.
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
    inhomogeneity should be printed in the terminal, at each step of the
    Add/Remove procedure.
  }
}
\details{
  All details of the Add/Remove procedure can be found in the article.
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
  \code{\link{BG.fitAll}}, \code{\link{BG.AddRemove}}
}
\examples{
## This example reproduces Figure 5 of the paper
## (Joucla et al., Cell Calcium, in revision)

## Create a single data set
## (fluorescence parameters are taken as the default arguments values)
set.seed(16051981)
pix.wrong <- c(5,10,15,20)
dB.pc <- 0.30
adu.list <- BG.createData(pix.wrong=pix.wrong, dB.pc=dB.pc, nRep=1)
adu <- adu.list[[1]]
tc <- attr(adu, "tc")
tOn <- attr(adu, "tOn")

# Fit a non-homogeneous-background model
directFit <- BG.track(adu, f.type = "mono", phi.type = "all",
                      pix.unif=NULL, tc, tOn, p_thresh=0.025,
                      Plot = "x11", Print = TRUE)
}
