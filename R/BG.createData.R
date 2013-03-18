BG.createData <-
function(nPixels=20,
                          phi=seq(0, 1, length.out=20+1)[-1],
                          Factor=500,
                          F0=1,
                          dF=0.5,
                          tau=1,
                          tOn=1,
                          tc=seq(0, 5, .1),
                          b=250,
                          dB.pc=0,
                          pix.wrong=NULL,
                          G=0.146,
                          nRep=1) {
  ## DATA (without noise)
  ## --------------------
  new.adu <- TRUE
  
  ## Choice of the pixel-specific scaling factors phi_i,
  ## as a round object (but square Region Of Interest)
  ## nPixels <- 20 ## ATTENTION ##
  ## phi <- seq(0, 1, length.out=nPixels+1)[-1]
  ## phi <- seq(0, 1, length.out=nPixels)
  phi <- round(phi,3)
  
  ## Global scaling factor
  ##Factor <- 500
  phi <- Factor * phi
  phi.backup <- phi
  
  ## Choice of the fluorescence transient f(t), normalize it
  ##F0 <- 1
  ##dF <- 0.5
  ##tau <- 1
  ##tOn <- 1
  ##t <- seq(0, 5, .1) ## ATTENTION ##
  f <- F0 + dF * ifelse(tc >= tOn, exp(-(tc-tOn)/tau), 0)
  
  ## Background fluorescence (homogeneous + wrong pixels)
  ##b <- 250
  ##if(!exists("dB.pc")) dB.pc <- 0
  ##if(!exists("pix.wrong")) pix.wrong <- NULL
  b.old <- b
  b.vec <- rep(b, nPixels)
  b.vec[pix.wrong] <- b.old * (1+dB.pc)
  B <- matrix(b.vec, nPixels, length(tc), byrow=!TRUE)
  
  ## Construction of the whole ideal adu matrix
  ## (including background fluorescence)
  adu <- phi %o% f + B
  adu.ideal <- adu
  
  ## True values of the parameters
  ## -----------------------------
  ## For f(t)
  if(!exists("f.type")) f.type <- "mono"
  if(f.type == "mono") {
    param.true <- c(dF = dF, tau = tau)
  } else {
    param.true <- c()
    for(fp in 2:length(f)) {
      S <- sprintf("param.true <- c(param.true, f_%d = f[%d])", fp, fp)
      eval(parse(text=S))
    }
  }
  ## For phi
  for(ph in 1:length(phi)) {
    S <- sprintf("param.true <- c(param.true, phi_%d = phi[%d])", ph, ph)
    eval(parse(text=S))
  }
  ## For b
  S <- sprintf("param.true <- c(param.true, b = %f)", b)
  eval(parse(text=S))
  
  ## For the b_i
  if(length(pix.wrong) > 0) {
    for(bi in pix.wrong) {
      S <- sprintf("param.true <- c(param.true, b_%d = B[%d,1])", bi, bi)
      eval(parse(text=S))
    }
  }
  
  ## Generate data according to the current noise model
  ## and the number of replicated data sets
  ## --------------------------------------------------
  ## G <- 0.146
  adu.list <- lapply(1:nRep,
                     function(i) {
                       adu.noise <- G * rpois(length(adu.ideal), adu.ideal / G)
                       adu <- matrix(adu.noise, nrow(adu.ideal), ncol(adu.ideal))
                       dimnames(adu) <- list(1:nrow(adu), NULL)
                       
                       if(new.adu == TRUE) {
                         ## Add the attributes to the adu
                         ## -----------------------------
                         attr(adu, "param.true") <- param.true
                         attr(adu, "f") <- f
                         attr(adu, "tc") <- tc
                         attr(adu, "tOn") <- tOn
                         attr(adu, "nPixels") <- nPixels
                       }                         
                       return(adu)
                     })
  
  ## adu <- adu + matrix(rnorm(length(adu), sd=50), nrow(adu), ncol(adu))
  
  ## adu.backup <- adu
  ## tc.backup <- tc
  
  ## -----------------------------
  return(adu.list)
}

