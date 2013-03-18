\name{POMC}
\alias{POMC}
\docType{data}
\title{
Calcium imaging data set obtained from POMC neurons of the mouse
hypothalamus
}
\description{
Experiments were performed on brain slices from 2-week-old male C57BL/6N
mice that expressed green fluorescent protein (eGFP) selectively in
proopiomelanocortin (POMC) neurons. POMC neurons were recorded at room
temperature in the whole-cell patch-clamp configuration. Calcium
transients were elicited by stepping the voltage-clamped membrane
potential from $-80$ to $+10$ mV for 50 ms. Fluorescence images were
acquired at 340 nm every 150 ms (exposure time: 12 ms). The imaging
setup consisted of an Imago/SensiCam CCD camera with a 640x480 chip
(Till Photonics, Grafelfing, Germany) and a Polychromator IV (Till
Photonics) that was coupled via an optical fiber into the upright
microscope. Emitted fluorescence was detected through a 440 nm long-pass
filter (LP440). Data were acquired as 80x60 frames using 8x8 on-chip
binning. Images were recorded in analog-to-digital units (ADUs) and
stored and analyzed as 12 bit grayscale images. The characteristics of
the CCD camera (gain, read-out noise) have been evaluated in a previous
study (Joucla et al., JNeurophysiol, 2010): G=0.146, standard deviation
of the read-out noise=16.4 photo-electrons.
}
\usage{data(POMC)}
\format{
  The format is:
List of 4
 $ adu.img: int [1:60, 1:80, 1:168] 285 285 294 292 291 291 279 298 290 294 ...
 $ time   : num [1:168] 525 525 525 526 526 ...
 $ tOn    : num 527
 $ ROI    :List of 2
  ..$ i: num [1:52] 25 25 25 25 26 26 26 26 26 27 ...
  ..$ j: num [1:52] 36 37 38 39 36 37 38 39 40 35 ...
}
\details{
Experiments were performed by Andreas Pippow, from the Kloppenburg lab
(University of Cologne, Germany).
}
\references{
Joucla S, Franconville R, Pippow A, Koppenburg P, Pouzat C (2013)
Estimating background-subtracted fluorescence transients in calcium
imaging experiments: a quantitative approach (Cell Calcium, in
revision)
}
\examples{
data(POMC)
str(POMC)
}
\keyword{datasets}
