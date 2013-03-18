BG.AddRemove <-
function(adu = NULL,
                         phi.type = "all",
                         f.type = "non-param",
                         tc = NULL,
                         tOn = NULL,
                         pix.unif = NULL,
                         pix.nonUnif = NULL,
                         nPixels = NULL,
                         p_thresh = 0.025,
                         Plot=FALSE,
                         Print=FALSE) {
    ## Track wrong pixels - subfunction for one step "Add/Remove"

    ## Initialize the stopProcedure variable
    stopProcedure <- FALSE
  
    ## --------
    ## Step ADD
    ## --------
    
    ## Initialise the Results variables
    bResults.val <- matrix(NA, nPixels, 2)
    bResults.icl <- matrix(NA, nPixels, 2)
    bResults.ich <- matrix(NA, nPixels, 2)
    bResults.pval <- rep(-0.02, nPixels)
    bResults.pc <- rep(NA, nPixels)

    ## Perform nPixels fits with only one single neuron with a different background
    for(k in pix.unif) {
        if(Print) print(sprintf("Add step: Test the background homogeneity of pixel %d/%d",k,length(pix.unif)))
        pix2test <- k
        pix.b.unif <- setdiff(pix.unif, k)
        
        ## Fit the data
        directFit_nh <- BG.fitAll(adu, phi.type, b.type="non-unif", pix.b.unif, f.type, tc, tOn)
        
        ## Update the Results variables
        bResults <- summary(directFit_nh)$parameters[c("b", sprintf("b_%d",k)),]
        bResults.val[k,] <- bResults[,1]
        bResults.icl[k,] <- bResults[,1] - bResults[,2]
        bResults.ich[k,] <- bResults[,1] + bResults[,2]
        bResults.pval[k] <- pnorm(-abs(bResults[2,1]-bResults[1,1]) /
                                  sqrt(bResults[1,2]^2+bResults[2,2]^2))
        if(Print) print(sprintf("  >>  its p-value is: %0.3f", bResults.pval[k]))
        bResults.pc[k] <- abs((bResults[2,1]-bResults[1,1])/bResults[1,2])
    }

    ## Plot all estimated b and p-values
    if(Plot == "pdf" | Plot == "x11" | Plot == "jpeg" | Plot == "tiff") {
        ## Use the current figure with the adequate layout
        
        ## Subplot1: estimated b and b_i for each pixel i
        plot(bResults.val[,1], pch=21, col="grey", bg="grey", xlim=c(0.5, nPixels+0.5),
             xaxp=c(-10, nPixels+10, 1), xlab="", xaxs="i",
             ## ylim=range(cbind(bResults.icl,bResults.ich), na.rm=TRUE),
             ylim=c(190,355), cex=0.6, mgp=c(3,0.75,0),
             ylab="", ## expression(paste("Estimated b and ", b[i])),
             yaxs="r", bty="l")
             ## main=sprintf("Step %d.1", count))
        lines(rep(1:nrow(bResults.val), each=3),
              as.vector(rbind(bResults.icl[,1], bResults.ich[,1],
                              rep(NA, nrow(bResults.ich)))), col="grey", bg="grey")
        points(bResults.val[,2], pch=21, col="black", bg="black", cex=0.6)
        lines(rep(1:nrow(bResults.val), each=3),
              as.vector(rbind(bResults.icl[,2], bResults.ich[,2],
                              rep(NA, nrow(bResults.ich)))), col="black")
        axis(side=1, at=c(0, 1, seq(0, nPixels-1, 5), nPixels, nPixels+1), mgp=c(3,0.6,0))

        ## Subplot2: estimated p-values
        plot(bResults.pval, pch=21, bg="white", col="white", xlim=c(0.5, nPixels+0.5),
             xaxs="i", yaxs="i", xaxp=c(-10, nPixels+10, 1), xlab="",
             ylim=c(-0.05,0.5), mgp=c(3,0.75,0), ylab="", ## "Estimated p-value",
             yaxp=c(0, 0.5, 5), bty="l")
             ## main=sprintf("Step %d.2", count))
        abline(h = c(0, p_thresh), lty=2, col=c("black", "red"))
        polygon(c(0.5, rep(nPixels+0.5,2), 0.5), rep(c(-0.05,0), each=2),
                col="grey", border="grey")
        points(pix.unif, bResults.pval[pix.unif], pch=21, col="black", bg="black", cex=0.6)
        points(pix.nonUnif, rep(-0.02, length(pix.nonUnif)), pch=25,
               col="red", bg="red", cex=0.6)
        axis(side=1, at=c(0, 1, seq(0, nPixels-1, 5), nPixels, nPixels+1), mgp=c(3,0.6,0))
    }

    ## Select the pixels that have a p-value below p_thresh,
    if(Print) print(bResults.pval)
    if(Print) print(p_thresh)
    pix2exclude <- which((bResults.pval < p_thresh & bResults.pval > 0)) ## |
                         ## (bResults.pval < 0.100 & bResults.pval > 0 & bResults.pc > 4))
    pix2exclude_old <- pix2exclude
    if(length(pix2exclude)==0) stopProcedure <- TRUE

    ## -----------
    ## Step REMOVE
    ## -----------
    
    ## Perform a direct fit to calculate more accurate b and p-values
    ## and control which ones are actually below p_thresh.
    print("Remove step: Test the background homogeneity of pixels: ")
    pix2exclude <- unique(sort(c(pix.nonUnif, pix2exclude)))
    print(pix2exclude)
    directFit_nh <- BG.fitAll(adu, phi.type, b.type="non-unif", pix.unif=setdiff(pix.unif, pix2exclude), f.type, tc, tOn)
    
    ## Update the Results variables
    bResults <- summary(directFit_nh)$parameters
    idx.b <- which(dimnames(bResults)[[1]] == "b")
    bResults <- bResults[-(1:(idx.b-1)),]
    bResults.val <- matrix(NA, nPixels, 2)
    bResults.icl <- matrix(NA, nPixels, 2)
    bResults.ich <- matrix(NA, nPixels, 2)
    bResults.pval <- rep(NA, nPixels)
    bResults.pc <- rep(NA, nPixels)

    for(p2e in 1:length(pix2exclude)) {
        bResults.val[pix2exclude[p2e],] <- bResults[c(1,p2e+1),1]
        bResults.icl[pix2exclude[p2e],] <- bResults[c(1,p2e+1),1] - bResults[c(1,p2e+1),2]
        bResults.ich[pix2exclude[p2e],] <- bResults[c(1,p2e+1),1] + bResults[c(1,p2e+1),2]
        bResults.pval[pix2exclude[p2e]] <- pnorm(-abs(bResults[p2e+1,1]-bResults[1,1]) /
                                                 sqrt(bResults[p2e+1,2]^2+bResults[1,2]^2))
        bResults.pc[pix2exclude[p2e]] <- abs((bResults[p2e+1,1]-bResults[1,1])/bResults[1,2])
    }

    ## Plot all estimated b and p-values (again...)
    if(Plot == "pdf" | Plot == "x11" | Plot == "jpeg" | Plot == "tiff") {
        ## Subplot3: estimated b and b_i for each pixel i in pix2exclude
        plot(bResults.val[,1], pch=21, col="white", bg="white",
             xlim=c(0.5, nPixels+0.5),
             xaxp=c(-10, nPixels+10, 1), xlab="azerty", xaxs="i",
             ## ylim=range(cbind(bResults.icl,bResults.ich), na.rm=TRUE),
             ylim=c(190,355),
             cex=0.6, mgp=c(3,0.75,0),
             ylab="", ## expression(paste("Estimated b and ", b[i])),
             yaxs="r", bty="l")
             ## main=sprintf("Step %d.3", count))
        abline(h = c(sort(unique(bResults.icl[,1])),
                     sort(unique(bResults.val[,1])),
                     sort(unique(bResults.ich[,1]))),
               lty = c(2,1,2), col="grey")
        points(bResults.val[,2], pch=21, col="black", bg="black", cex=0.6)
        lines(rep(1:nrow(bResults.val), each=3),
              as.vector(rbind(bResults.icl[,2], bResults.ich[,2],
                              rep(NA, nrow(bResults.ich)))), col="black")
        axis(side=1, at=c(1, seq(0, nPixels-1, 5), nPixels), mgp=c(3,0.6,0))
        ## Subplot4: estimated p-values
        plot(bResults.pval, pch=21, bg="white", col="white",
             xlim=c(0.5, nPixels+0.5),
             xaxs="i", yaxs="i", xaxp=c(-10, nPixels+10, 1), xlab="",
             ylim=c(-0.05,0.5), mgp=c(3,0.75,0),
             ylab="", ## "Estimated p-value",
             yaxp=c(0, 0.5, 5), bty="l")
             ## main=sprintf("Step %d.4", count))
        abline(h = c(0, p_thresh), lty=2, col=c("black", "red"))
        points(bResults.pval, pch=21, col="black", bg="black", cex=0.6)
        axis(side=1, at=c(0, 1, seq(0, nPixels-1, 5), nPixels, nPixels+1), mgp=c(3,0.6,0))
    }

    ## Update pix2exclude
    if(length(pix2exclude) >= 1) {
      pix2exclude <- which((bResults.pval < p_thresh & bResults.pval > 0))
    }
    if(length(intersect(pix2exclude_old, pix2exclude)) == 0) {
        pix2exclude <- NULL
    }

    return(pix2exclude)
}

