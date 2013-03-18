BG.fluoNormFct <-
function(f.type = "non-param",
                           nTimePoints = NULL
                           ) {
  ## Create an empty function
  fluoNormFct <- function() NULL
  
  ## Define its arguments ...
  S4F <- "formals(fluoNormFct) <- list("
  if(f.type == "non-param") {
    f_1 <- 1 ## The first fluorescence value is set to 1
    idx.f <- 2:(nTimePoints-1)
    for(j in idx.f) {
      S4F <- c(S4F, sprintf("f_%d=NULL,",j))
    }
    S4F <- c(S4F, sprintf("f_%d=NULL)",nTimePoints))
  } else if(f.type == "mono") {
    S4F <- c(S4F, "tc=NULL,")
    S4F <- c(S4F, "tOn=NULL,")
    S4F <- c(S4F, "dF=NULL,")
    S4F <- c(S4F, "tau=NULL)")
  } else if(f.type == "fct_b") {
    ## Close the parenthesis
    S4F <- c(S4F, ")")
  }
  eval(parse(text=S4F))
  
  ## ... and its body
  S4B <- "body(fluoNormFct) <- expression({"
  if(f.type == "non-param") {
    ## S4B <- c(S4B, "f <- 1;")
    ## for(j in idx.f) {
    ##   S4B <- c(S4B, sprintf("f <- c(f, f_%d);",j))
    ## }
    ## S4B <- c(S4B, sprintf("f <- c(f, f_%d);",nTimePoints))
    S4B <- c(S4B, "f <- c(1,")
    for(j in idx.f) {
      S4B <- c(S4B, sprintf("f_%d,",j))
    }
    S4B <- c(S4B, sprintf("f_%d)", nTimePoints))
  } else if(f.type == "mono") {
    S4B <- c(S4B, sprintf("f <- 1 + dF * ifelse(tc >= tOn, exp(-(tc-tOn)/tau), 0);"))
  } else if(f.type == "fct_b") {
    ## Do nothing
  }
  S4B <- c(S4B, "})")
  eval(parse(text=S4B))

  return(fluoNormFct)
}

