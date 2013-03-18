BG.fitAll <-
function(adu = NULL,
                      phi.type = "all",
                      b.type="unif",
                      pix.unif = NULL,
                      f.type = NULL,
                      tc = NULL,
                      tOn = NULL
                      ) {
  ## Initial guesses
  ig <- BG.ig(adu, phi.type, f.type, b.type, pix.unif)
  
  ## Fit function
  SQRT <- TRUE
  fitFct <- BG.mkFct(adu, phi.type, b.type, pix.unif, f.type, SQRT)

  ## Define the function for the fit
  Formals <- names(formals(fitFct))
  myfunction <- "fitFct("
  for(i in Formals[-length(Formals)]) myfunction <- c(myfunction, i, ",")
  myfunction <- c(myfunction, Formals[length(Formals)], ")")
  
  ## Perform the fit
  adu4Fit <- as.vector(adu^(1/(1+SQRT)))
  directFit <- c()
  S4Fit <- c("directFit <- nls(adu4Fit ~ ", myfunction)
  S4Fit <- c(S4Fit,",start = ig,")
#   S4Fit <- c(S4Fit,",trace = TRUE,")
  S4Fit <- c(S4Fit,"control = list(warnOnly=!TRUE))")
  eval(parse(text=S4Fit))
  
  return(directFit)
}

