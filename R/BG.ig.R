BG.ig <-
function(adu=NULL,
                  phi.type="all",
                  f.type="non-param",
                  b.type="unif",
                  pix.unif=NULL
                  ) {
  ## We will determine initial guesses for all phi_i, all f_j and b
  ## based on the whole adu matrix
  adu.mean.time <- rowMeans(adu)
  adu.mean.space <- colMeans(adu)
  adu.mean <- mean(adu)

  ## phi.mean represents the mean of all phi_i
  bOfPhi.Mean <- function(phi.mean) adu.mean - phi.mean
  fOfPhi.Mean <- function(phi.mean) (adu.mean.space - adu.mean) / phi.mean + 1
  ## Here, we assume that the time course of f(t) has a normalized integral (=1)
  phiOfPhi.Mean <- function(phi.mean) adu.mean.time - adu.mean + phi.mean

  ## Find an initial guess for phi.mean
  if(!exists("b.drawn")) b.drawn <- 0
  if(!exists("tc")) tc <- seq(0,5,0.1)
  if(b.type == "unif" | b.type == "non-unif") {
    LOFPhi.Mean <- function(phi.mean) {
      sum((adu - phiOfPhi.Mean(phi.mean) %o% fOfPhi.Mean(phi.mean) - bOfPhi.Mean(phi.mean))^2)
    }
  } else if(b.type == "none" | b.type == "one") {
    LOFPhi.Mean <- function(phi.mean) {
      sum((adu - phiOfPhi.Mean(phi.mean) %o% fOfPhi.Mean(phi.mean) - mean(b.drawn))^2)
    }
  }
  
  ## Fit the value of phi.mean
  phi.mean <- optimize(LOFPhi.Mean, interval=c(0, sum(adu)))$minimum
  
  ## Initial guesses (don't forget to "unnormalize" f and phi...
  if(b.type == "unif" | b.type == "non-unif") {
    b <- bOfPhi.Mean(phi.mean)
  } else if(b.type == "none") {
    b <- mean(b.drawn)
  }
  f <- fOfPhi.Mean(phi.mean)
  phi <- phiOfPhi.Mean(phi.mean)
  phi <- phi * f[1]
  f <- f / f[1]
  phi[phi<0] <- 0.001 ## All phi_i must be positive

  ## Create the string for the creation of the ig variable
  ig <- c()
  if(!exists("tOn")) tOn <- 1
  S4IG <- "ig <- list("
  
  ## 1) The f(t) time course
  if(f.type == "non-param") {
    f_1 <- 1 ## The first fluorescence value is set to 1
    idx.f <- 2:ncol(adu)
    for(j in idx.f) {
      eval(parse(text=sprintf("f_%d <- f[%d]",j,j)))
      S4IG <- c(S4IG, sprintf("f_%d = f_%d,",j,j))
    }
  } else if(f.type == "mono") {
    ## Create the function that should fit the fluorescence time course
    fluoNormFct <- BG.fluoNormFct(f.type, ncol(adu))
    
    ## Fit the time course
    f.cut <- f[which(tc>tOn)]/f[1]-1
    tc.cut <- (tc-min(tc))[which(tc>tOn)]
    tc.cut <- tc.cut - min(tc.cut)
    L <- lm(log(f.cut) ~ tc.cut)
    if(tc.lim <- -3/coef(L)[2] < max(tc.cut)) {
      f.cut <- f.cut[1:which(tc.cut >= tc.lim)[1]]
      tc.cut <- tc.cut[1:which(tc.cut >= tc.lim)[1]]
      L <- lm(log(f.cut) ~ tc.cut)
    }
    dF.est <- exp(coef(L)[1])
    tau.est <- -1/coef(L)[2]
    
    S4IG <- c(S4IG, sprintf("dF = %f,", dF.est))
    S4IG <- c(S4IG, sprintf("tau = %f,", tau.est))
  } else if(f.type == "fct_b") {
    ## Do nothing
  }
  
  ## 2) The phi_i(s) pixel-specific amplitude factors
  idx.phi <- 1:nrow(adu)
  if(phi.type == "all") {
    for(i in idx.phi[-length(idx.phi)]) {
      eval(parse(text=sprintf("phi_%d <- phi[%d]",i,i)))
      S4IG <- c(S4IG, sprintf("phi_%d = phi_%d,",i,i))
    }
    i <- length(idx.phi)
    eval(parse(text=sprintf("phi_%d <- phi[%d]",i,i)))
    if(b.type != "none") {
      S4IG <- c(S4IG, sprintf("phi_%d = phi_%d,",i,i))
    } else if(b.type == "none") {
      S4IG <- c(S4IG, sprintf("phi_%d = phi_%d)",i,i))
    }
  } else if(phi.type == "fct_b") {
    ## Do nothing
  }
  
  ## 3) The background fluorescence value(s) b (and b_i(s))
  if(b.type == "non-unif" & length(setdiff(idx.phi, pix.unif)) >= 1) {
    idx.b <- 1:nrow(adu)
    if(length(pix.unif) > 1) idx.b <- idx.b[-pix.unif]
    S4IG <- c(S4IG, "b = b,")
    for(i in idx.b[-length(idx.b)]) {
      eval(parse(text=sprintf("b_%d <- b;",i)))
      S4IG <- c(S4IG, sprintf("b_%d = b_%d,",i,i))
    }
    eval(parse(text=sprintf("b_%d <- b;",idx.b[length(idx.b)])))
    S4IG <- c(S4IG, sprintf("b_%d = b_%d)",idx.b[length(idx.b)],idx.b[length(idx.b)]))
  } else {
    S4IG <- c(S4IG, "b = b)")
  }
  
  ## Create the variable ig and return it
  eval(parse(text=S4IG))
  return(ig)
}

