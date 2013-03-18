BG.mkFct <-
function(adu=NULL,
                     phi.type="all",
                     b.type="unif",
                     pix.unif=NULL,
                     f.type = NULL,
                     SQRT=TRUE
                     ) {
  ## Define the fluoNormFct function
  fluoNormFct <- BG.fluoNormFct(f.type, ncol(adu))
  
  ## Create an empty function
  fitFct <- function() NULL
  
  ## ------------------------
  ## Define its arguments ...
  ## ------------------------
  S4F <- "formals(fitFct) <- list("
  
  ## 1) for the f(t) normalized fluorescence time course
  Formals <- names(formals(fluoNormFct))
  for(j in Formals) {
    S4F <- c(S4F, sprintf("%s=NULL,", j))
  }
  
  ## 2) for the phi_i(s) pixel-specific amplitude factors
  idx.phi <- 1:nrow(adu)
  if(phi.type == "all") {
    for(i in idx.phi[-length(idx.phi)]) {
      S4F <- c(S4F, sprintf("phi_%d=NULL,",i))
    }
    i <- length(idx.phi)
    if(b.type != "none") {
      S4F <- c(S4F, sprintf("phi_%d=NULL,",i))
    } else {
      ## Stop the list of arguments after the last phi_i
      S4F <- c(S4F, sprintf("phi_%d=NULL)",i))
    }
  } else if(phi.type == "fct_b") {
    ## Do nothing
  }
  
  ## 3) for the background fluorescence value(s) b (and b_i(s))
  if(b.type == "non-unif" & length(setdiff(idx.phi, pix.unif)) >= 1) {
    idx.b <- 1:nrow(adu)
    if(length(pix.unif) > 1) {
      idx.b <- idx.b[-pix.unif]
    }
    S4F <- c(S4F, "b=NULL,")
    for(i in idx.b[-length(idx.b)]) {
      S4F <- c(S4F, sprintf("b_%d=NULL,",i))
    }
    S4F <- c(S4F, sprintf("b_%d=NULL)",idx.b[length(idx.b)]))
  } else if(b.type == "unif") {
    S4F <- c(S4F, "b=NULL)")
  }else if(b.type == "none") {
    ## the b_i(s) are not parameters anymore, but should be variables in the environment
  }

  eval(parse(text=S4F))

  ## ----------------  
  ## ... and its body
  ## ----------------  
  S4B <- "body(fitFct) <- expression({"

  ## 1) for the background
  if(b.type == "non-unif" & length(setdiff(idx.phi, pix.unif)) >= 1) {
    ## Define the background as the concatenation of b and the b_i(s)
    idx.b <- 1:nrow(adu)
    S4B <- c(S4B, "b <- c(")
    for(i in idx.b[-length(idx.b)]) {
      if(length(which(pix.unif == i)) == 1) {
        S4B <- c(S4B, "b,")
      } else {
        S4B <- c(S4B, sprintf("b_%d,",i))
      }
    }
    i <- idx.b[length(idx.b)]
    if(length(which(pix.unif == i)) == 1) {
      S4B <- c(S4B, "b);")
    } else {
      S4B <- c(S4B, sprintf("b_%d);",i))
    }
    N <- ncol(adu)
    S4B <- c(S4B, sprintf("b <- matrix(rep(b, %d), ncol=%d);", N, N))
  } else if(b.type == "unif") {
    ## Do nothing
    ## S4B <- c(S4B, "b <- matrix(b, nrow(adu), ncol(adu));")
  } else if(b.type == "none") {
    ## Create a "b" vector with all b_i values concatenated
    ## (the b_i values are in the b.drawn vector)
    idx.b <- 1:nrow(adu)
    S4B <- c(S4B, "b <- b.drawn;")
    N <- ncol(adu)
    S4B <- c(S4B, sprintf("b <- matrix(rep(b, %d), ncol=%d);", N, N))
  }
  
  ## 2) for the f(t) fluorescence time course
  Formals <- names(formals(fluoNormFct))
  if(length(Formals) > 0) {
    S4B <- c(S4B, "f <- fluoNormFct(")
    for(j in Formals[-length(Formals)]) {
      S4B <- c(S4B, sprintf("%s,", j))
    }
    S4B <- c(S4B, sprintf("%s);", Formals[length(Formals)]))
  } else { ## This case corresponds to the case f.type == "fct_b"
    ## Calculate the average of all background fluorescence values
    if(b.type == "unif") {
      S4B <- c(S4B, "f <- (colMeans(adu)-b)/(mean(adu[,1])-b);")
    } else {
      S4B <- c(S4B, "b.mean <- mean(b);")
      S4B <- c(S4B, "f <- (colMeans(adu)-b.mean)/(mean(adu[,1])-b.mean);")
    }
  }
  
  ## for the phi_i's then
  if(phi.type == "all") {
    S4B <- c(S4B, "phi <- c();")
    for(i in idx.phi) {
      S4B <- c(S4B, sprintf("phi <- c(phi, phi_%d);",i))
    }
  } else if(phi.type == "fct_b") {
    if(b.type == "unif") {
      S4B <- c(S4B, "phi <- (rowMeans(adu) - b) / mean(f);")
    } else if(b.type == "none") {
      S4B <- c(S4B, "phi <- c();")
      for(k in idx.phi) {
        S4B <- c(S4B, sprintf("phi <- c(phi, (mean(adu[%d,]) - b_%d) / mean(f));", k, k))
      }
    } else if(b.type == "non-unif") {
      S4B <- c(S4B, sprintf("phi <- vector('numeric', %d);", length(idx.phi)))
      for(k in idx.phi)  {
        if(length((which(idx.b == k))) == 1) {
          S4B <- c(S4B, sprintf("phi[%d] <- (mean(adu[%d,]) - b_%d) / mean(f);", k, k, k))
        } else {
          S4B <- c(S4B, sprintf("phi[%d] <- (mean(adu[%d,]) - b) / mean(f);", k, k))
        }
      }
    }
  }
  
  ## Create the whole adu matrix
  if(SQRT == TRUE) {
    S4B <- c(S4B, "adu <- sqrt(as.vector(phi %o% f + b));")
  } else {
    S4B <- c(S4B, "adu <- as.vector(phi %o% f + b);")
  }
  
  S4B <- c(S4B, "return(adu)})")
  eval(parse(text=S4B))

  return(fitFct)
}

