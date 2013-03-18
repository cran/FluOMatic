\name{BG.fluoNormFct}
\alias{BG.fluoNormFct}
\title{
Create a function used to calculate a normalized fluorescence transient.
}
\description{
This low-level function should be used to create a function calculating
a normalized fluorescence transient. The latter is used in other
functions of the package, such as \code{BG.ig} or \code{BG.mkFct}.
}
\usage{
BG.fluoNormFct(f.type = "non-param", nTimePoints = NULL)
}
\arguments{
  \item{f.type}{
    A string specifying if the normalized fluorescence transient f(t) is
    modeled by a monoexponential model (set \code{f.type} to "mono") or
    is estimated pointwise (set \code{f.type} to "non-param").
  }
  \item{nTimePoints}{
    A numerical value giving the length (in time samples) of the desired
    noramlized fluorescence transient (should be specified only if
    \code{f.type} is set to "non-param").
  }
}
\value{
  A function with arguments:
  \item{tc, tOn, dF, tau }{if \code{f.type} is set to "mono". In that
    case, the body of the output function creates a vector with a
    monoexponential time course.}
  \item{f_2, f_3, ..., f_N}{if \code{f.type} is set "non-param" (N
    corresponds to \code{nTimePoints}). In that case, the body of the
    output function creates a vector with first value equal to 1 (in our
    paradigm, the normalized fluorescence baseline is equal to 1),
    followed by all inuput arguments.}
}
\author{
  Sebastien Joucla \email{joucla_sebastien@yahoo.fr}
}
\seealso{
  \code{\link{BG.ig}}, \code{\link{BG.mkFct}}
}
\examples{
## Basic example creating a function corresponding either to the
## monoexponential case or to the non-parametric case
F1 <- BG.fluoNormFct(f.type = "mono")
F2 <- BG.fluoNormFct(f.type = "non-param", nTimePoints=51)
}
