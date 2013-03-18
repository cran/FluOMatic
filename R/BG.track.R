BG.track <- function(adu = NULL,
                     f.type = "mono",
                     phi.type = "all",
                     pix.unif = NULL,
                     tc = NULL,
                     tOn = NULL,
                     p_thresh = 0.025,
                     Plot = "x11",
                     Print = FALSE
                     ) {
    ## Function BG.track
    ## 
    ## This function fits 1 set of data with non-uniform background
    ## with a model depending on f.type and phi.type.
    ## This is a recursive approach aiming at finding, at each step,
    ## a supplementary pixel with a background different than the others.
    ## The function return the last "fit" object
    ## with the number of "different" pixels as attribute.
    
    ## Fit the whole model on all data (uniform background)
    ## fluoNormFct <- BgNonUnif_fluoNormFct(f.type, ncol(adu))
    b.type <- "unif"
    directFit <- BG.fitAll(adu, phi.type, b.type, pix.unif, f.type, tc, tOn)
    attr(directFit, "RawData") <- adu
    class(directFit) <- c("BgUnifFitObject", "nls")
    nParam <- nrow(summary(directFit)$parameters)
    nParam_butB <- nParam-1

    ## Plot the raw and predicted data for each pixel
    if(FALSE) plot(directFit)

    ## ----------------------------------------
    ## Compute b and b_i for pixels i that are
    ## suspected to have a different background
    ## ----------------------------------------

    ## Step1: Start with all pixels homogeneous except 1
    ## compute a p-value for each of them
    nPixels <- nrow(adu)
    if(length(pix.unif) == 0) pix.unif <- 1:nPixels
    pix.unif.bk <- pix.unif
    pix.nonUnif <- c()
    count <- 1
        
    ## New: 27/08/2012
    ## Prepare the figure for the plots of each add/remove steps
    if(Plot == "pdf" | Plot == "x11" | Plot == "jpeg" | Plot == "tiff") {
        if(Plot == "pdf") {
          pdf("Figure5.pdf", width=7, height=4.5)
        } else if(Plot == "jpeg") {
          jpeg("Figure5.jpg", width=7, height=4.5, units="in", res=300, quality=90)          
        } else if(Plot == "tiff") {
          tiff("Figure5.tiff", width=7, height=4.5, units="in", res=300, compression="lzw")          
        } else {
          x11(width=7, height=4.5)
        }
        layout(rbind(cbind(1,2,2,3,3), cbind(4,5,6,7,8), cbind(9,10,11,12,13), cbind(14,15,16,17,18), cbind(19,20,21,22,23)),
               heights=c(0.05,0.035,0.305,0.305,0.305)*7, widths=c(0.035,0.23125,0.23125,0.23125,0.23125))
        
        par(mar=rep(0,4), oma=rep(0,4))
        par(las=1)
    
        ##Top Left Corner: nothing
        plot(NA, col='white', ann=FALSE, axes=FALSE, ylim=c(0,1), xlim=c(0,1))
        
        ## Top Middle: "add"
        plot(NA, col='white', ann=FALSE, axes=FALSE, ylim=c(0,1), xlim=c(0,1), xaxs="i", yaxs="i")
        polygon(c(0,1,1,0),c(0,0,1,1), col="gray", border="white", lwd=4, ylim=c(0,1))
        text(0.5, 0.5, '"Add"', col="black", font=2, cex=1.2)
        
        ## Top Right: "remove"
        plot(NA, col='white', ann=FALSE, axes=FALSE, ylim=c(0,1), xlim=c(0,1), xaxs="i", yaxs="i")
        polygon(c(0,1,1,0),c(0,0,1,1), col="gray", border="white", lwd=4, ylim=c(0,1))
        text(0.5, 0.5, '"Remove"', col="black", font=2, cex=1.2)
        
        ## Main titles
        plot(NA, col='white', ann=FALSE, axes=FALSE, ylim=c(0,1), xlim=c(0,1), xaxs="i", yaxs="i")
        
        plot(NA, col='white', ann=FALSE, axes=FALSE, ylim=c(0,1), xlim=c(0,1), xaxs="i", yaxs="i")
        text(0.5, 0.5, expression(paste("  Estimated b and ", b[i])), col="black", cex=1.2)
        
        plot(NA, col='white', ann=FALSE, axes=FALSE, ylim=c(0,1), xlim=c(0,1), xaxs="i", yaxs="i")
        text(0.5, 0.5, expression(paste("     Estimated ", p[i], "-values")), col="black", cex=1.2)
        
        plot(NA, col='white', ann=FALSE, axes=FALSE, ylim=c(0,1), xlim=c(0,1), xaxs="i", yaxs="i")
        text(0.5, 0.5, expression(paste("  Estimated b and ", b[i])), col="black", cex=1.2)
        
        plot(NA, col='white', ann=FALSE, axes=FALSE, ylim=c(0,1), xlim=c(0,1), xaxs="i", yaxs="i")
        text(0.5, 0.5, expression(paste("     Estimated ", p[i], "-values")), col="black", cex=1.2)
        
        ## Line2, Left: Step 1
        plot(NA, col='white', ann=FALSE, axes=FALSE, ylim=c(0,1), xlim=c(0,1), xaxs="i", yaxs="i")
        polygon(c(0,1,1,0),c(0,0,1,1), col="gray", border="white", lwd=4, ylim=c(0,1))
        text(0.5, 0.5, 'Step 1', col="black", font=2, cex=1.2, srt=90)
    }
    
    ## --------------------
    ## Call the subfunction
    ## --------------------
    
    ## Create the fluoNormFunction
    ## fluoNormFct <- BG.fluoNormFct(f.type, ncol(adu))
    
    ## Step 1
    step_counter <- 1
    if(Plot == "pdf" | Plot == "x11" | Plot == "jpeg" | Plot == "tiff") par(mar=c(2,2.5,1,1))
    pix2exclude <- BG.AddRemove(adu, phi.type, f.type, tc, tOn,
                                pix.unif, pix.nonUnif, nPixels, p_thresh, Plot, Print)
        
    ## Following steps
    while(length(pix2exclude) != 0) {
        step_counter <- step_counter +1
      
        ## PLOT
        ## Line i+2, Left: Step i+1
        if(Plot == "pdf" | Plot == "x11" | Plot == "jpeg" | Plot == "tiff") {
            par(mar=rep(0,4))
            plot(NA, col='white', ann=FALSE, axes=FALSE, , ylim=c(0,1), xlim=c(0,1), xaxs="i", yaxs="i")        
            polygon(c(0,1,1,0),c(0,0,1,1), col="gray", border="white", lwd=4, ylim=c(0,1))
            text(0.5, 0.5, paste('Step ', step_counter), col="black", font=2, cex=1.2, srt=90)
        }
             
        ## Exclude them from pix.unif and add them to pix.nonUnif
        pix.unif <- setdiff(pix.unif.bk, pix2exclude)
        pix.nonUnif <- pix2exclude
        count <- count+1

        ## Step2: repeat step1 with the new pix.unif and pix.nonUnif vectors
        if(Plot == "pdf" | Plot == "x11" | Plot == "jpeg" | Plot == "tiff") par(mar=c(2,2.5,1,1))
        pix2exclude <- BG.AddRemove(adu, phi.type, f.type, tc, tOn,
                                    pix.unif, pix.nonUnif, nPixels, p_thresh, Plot, Print)
    }
    if(Plot == "pdf" | Plot == "jpeg" | Plot == "tiff") dev.off()
    
    ## Perform a last fit with only uniform pixels having a uniform background
    b.type <- "non-unif"
    directFit <- BG.fitAll(adu, phi.type, b.type, pix.unif, f.type, tc, tOn)
    
    attr(directFit, "pix.unif") <- pix.unif
    attr(directFit, "pix.nonUnif") <- pix.nonUnif
    ## attr(directFit, "adu.raw") <- adu
    class(directFit) <- c("nls", "BgUnifFitObject")

    return(directFit)
}
