
# Plot Permutation Test Results
#
# @description
#' Function for plotting the results from a \code{permTestResults} object.
#' 
#' @method plot permTestResults
#' 
#' @param x          an object of class \code{permTestResults}.
#' @param pvalthres  p-value threshold for significance. Default is 0.05.
#' @param plotType   the type of plot to display. This must be one of \code{"Area"} or \code{"Tailed"}. Default is \code{"Area"}.
#' @param main       a character specifying the title of the plot. Defaults to "".
#' @param xlab       a character specifying the label of the x axis. Defaults to NULL, which produces a plot with the evaluation function name as the x axis label.
#' @param ylab       a character specifying the label of the y axis. Defaults to "".
#' @param ylim       defines the y limits of the plot. Passed to the underlying \code{plot} call.
#' @param xlim       defines the x limits of the plot. Passed to the underlying \code{plot} call.
#' @param ...        further arguments to be passed to or from methods.
#' 
#' @return A plot is created on the current graphics device.
#' 
#' @seealso \code{\link{permTest}}
#' @import graphics
#' @importFrom stats dnorm qnorm rnorm runif
#'
#' @export

# hist(xcoords, prob = TRUE, ylim = ylim, breaks = 30, xlim = xlim,las = 1, col = "#d5e8eb", border = "#d5e8eb", xlab=xlab, ylab=ylab, main=paste(main, "p-value: ",pval, "  Z-score: ", zscore, "  n perm: ", nperm),cex.main=0.7,font.main=3, ...)


plot.permTestResults<-function(x, pvalthres=0.05, plotType="Tailed", main="", xlab=NULL, ylab="", ylim=NULL, xlim=NULL, ...){
  
  old.scipen <- options()$scipen
  
  options(scipen=999)
  
  if(class(x)!="permTestResults")  stop("x must be a permTestResults object")
  if(!is.numeric(pvalthres)) stop("pvalthres must be numeric")
  plotType<-match.arg(plotType,c("Area","Tailed"))
  
  
  if(is.null(xlab)) xlab <- paste0(x$evaluate.function.name)
  if(nchar(main)>0) main <- paste0(main, "\n")
  
  alternative<-x$alternative
  xcoords<-x$permuted
  xcoords<-xcoords[order(xcoords)]
  pval<-round(x$pval,4)
  nperm<-x$ntimes
  mperm<-mean(xcoords,na.rm=TRUE)
  mobs<-x$observed
  zscore<-round(x$zscore,3)
  
  if(is.finite(zscore)){
    y<-dnorm(xcoords,mean=mean(xcoords,na.rm=TRUE),sd=stats::sd(xcoords,na.rm=TRUE))
    xhist<-hist(xcoords,breaks=30,plot=FALSE)$density
    ymax<-max(max(y,na.rm=TRUE),max(xhist,na.rm=TRUE))
    
    if (alternative=="greater") aux<-qnorm((1-pvalthres),mean=mean(xcoords,na.rm=TRUE),sd=sd(xcoords,na.rm=TRUE))
    if (alternative=="less") aux<-qnorm(pvalthres,mean=mean(xcoords,na.rm=TRUE),sd=sd(xcoords,na.rm=TRUE))
    
    xmin<-min(mobs, min(xcoords,na.rm=TRUE), min(aux,na.rm=TRUE), na.rm=TRUE)
    xmax<-max(mobs, max(xcoords,na.rm=TRUE), max(aux,na.rm=TRUE), na.rm=TRUE)
    
    if(is.null(ylim)) ylim <- c(0,ymax)
    if(is.null(xlim)) xlim <- c(xmin,xmax)
    
    hist(xcoords, prob = TRUE, ylim = ylim, breaks = 30, xlim = xlim,
         las = 1, col = "#d5e8eb", border = "#d5e8eb", xlab=NULL, ylab=NULL, main=NULL,...)
    title(paste(main, "p-value: ",pval, "  Z-score: ", zscore, "  n perm: ", nperm), line =0.4, cex.main=1,font.main=3)
    title(xlab=xlab, line = 2)
    
    if(plotType=="Area"){
      if(alternative=="greater"){
        polygon(c(aux,aux,xmax,xmax),c(max(y,na.rm=TRUE),0,0,max(y,na.rm=TRUE)),col="#FD6D6D",density=10,border="white")
        lines(c(aux,aux),c(0,ymax*0.8),col="#FD6D6D",lwd=3)
        text(aux,ymax*0.9,bquote(alpha==.(pvalthres)),cex=0.8,pos=4)
      }
      if(alternative=="less"){
        polygon(c(aux,aux,xmin,xmin),c(max(y,na.rm=TRUE),0,0,max(y,na.rm=TRUE)),col="#FD6D6D",density=10,border="white")
        lines(c(aux,aux),c(0,ymax*0.8),col="#FD6D6D",lwd=3)
        text(aux,ymax*0.9,bquote(alpha==.(pvalthres)),cex=0.8)
      }
    }
    
    #I only modify tailed because I'm lazy, it is the default setting for plotPerm#
    if(plotType=="Tailed"){
      if(alternative=="greater"){
        aux3<-seq(aux,xmax,length=50)
        y3<-dnorm(aux3,mean(xcoords,na.rm=TRUE),sd(xcoords,na.rm=TRUE))
        polygon(c(aux3[1],aux3,aux3[length(aux3)]),c(0,y3,0),col="#FD6D6D",density=30,border="white")
        lines(c(aux,aux),c(0,ymax*0.8),col="#FD6D6D",lwd=3)
        #text(aux,ymax*0.3,bquote(alpha==.(pvalthres)),cex=0.8)
      }
      if(alternative=="less"){
        aux3<-seq(aux,xmin,length=50)
        y3<-dnorm(aux3,mean(xcoords,na.rm=TRUE),sd(xcoords,na.rm=TRUE))
        polygon(c(aux3[1],aux3,aux3[length(aux3)]),c(0,y3,0),col="#FD6D6D",density=30,border="white")
        lines(c(aux,aux),c(0,ymax*0.8),col="#FD6D6D",lwd=2)
        #text(aux,ymax*0.3,bquote(alpha==.(pvalthres)),cex=0.8)
      }
    }
    
    
    lines(xcoords,y,lwd=2)
    lines(c(mperm,mperm),c(0,ymax*0.8),col="black",lwd=3, lty=2)
    #text(mperm,ymax*0.9,expression(Ev[perm]),cex=0.8)
    
    lines(c(mobs,mobs),c(0,ymax*0.8),col="#38B3B7",lwd=3)
    text(mobs,ymax*0.9,expression(Ev[obs]),cex=1)
    arrows(mperm,ymax*0.75,mobs,ymax*0.75,length=0.1,code=3)
    box(lwd=1.2)
  }
  
  if(!is.finite(zscore)){
    xhist<-hist(xcoords,breaks=30,plot=FALSE, col = "#f58c8c")$density
    ymax<-max(xhist,na.rm=TRUE)
    
    if(is.null(ylim)) ylim <- c(0,ymax)
    if(is.null(xlim)) xlim <- c(min(mobs,min(xcoords,na.rm=TRUE),na.rm=TRUE), max(mobs,max(xcoords,na.rm=TRUE),na.rm=TRUE))
    
    #HISTOGRAM parameter here#
    hist(xcoords, prob = TRUE, ylim = ylim, breaks = 30, xlim = xlim, las = 1, col = "#f58c8c", border = "#f58c8c", xlab=xlab, ylab=ylab, main=paste(main, "p-value: ", pval, "\n Z-score: ", zscore, "\n n perm: ", nperm, "\n randomization: ", paste0(x$randomize.function.name)), cex.main=0.8, ...)
    
    lines(c(mperm,mperm),c(0,ymax*0.8),col="black",lwd=3)
    #text(mperm,ymax*0.9,expression(Ev[perm]),cex=0.8)
    
    lines(c(mobs,mobs),c(0,ymax*0.8),col="forestgreen",lwd=3)
    text(mobs,ymax*0.9,expression(Ev[obs]),cex=1)
    arrows(mperm,ymax*0.75,mobs,ymax*0.75,length=0.1,code=3)
    box(lwd=1.2)
    
    text(paste("p-value: ",pval, "  Z-score: ", zscore, "  n perm: ", nperm), ymax*1.1)
    
    warning(paste0("all permuted values are equal to ",xcoords[1],". It is not posible to adjust a normal distribution nor to compute a Z-score."))
    
  }
  
  options(scipen=old.scipen)
  
}


