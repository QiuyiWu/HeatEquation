## this file contains all the R functions required for a kernel
## smoother. Do *NOT* mix simulations with functions!!

## ks() is the main function implemented by Xing. Note that the main
## purpose for writing this function is to demonstrate the underlying
## mathematics, not for production-level use. For example, the current
## function is based on the very inefficient for-loop in R; we almost
## certainly will need to replace it by C/C++ in a production-level R
## package. See ksmooth.c included in R's source code for a reference
## implementation. Besides, it may be more efficient to implement KS
## based on fast discrete Fourier transformation for fixed bandwidth.
## Besides, for a production-level R package, we also need to add many
## QC code to it, for example, it should stop if length(x) !=
## length(y); it should have a mechanism to deal with NAs and probably
## outliers; it also needs to check whether x is ordered or not -- KS
## does not work for unsorted x!

## x is the *ordered* vector of the covariate, y is the observed value
## to be smoothed. bandwidth can either be a constant bandwidth or a
## vector (adaptive bandwidth) that has the same length as x. In both
## cases, the value (values) of bandwidth are defined as the standard
## deviation of the Gaussian kernel. min.bandwith and max.bandwidth
## are the minimum and maximum bandwidth. If bandwidth[i]
## <min.bandwidth, we simply return y[j] as yhat[j]; if bandwidth[i] >
## max.bandwidth, we use simple unweighted average over the smoothing
## window as yhat[j]. By default ("auto"), min.bandwidth is set to be
## min(xsteps)/4, where xsteps are the differences of two consecutive
## x[j]; max.bandwidth is set to be 10*xrange, where xrange is the
## range of the domain.

## max.window is the maximum smoothing window, by default, it is set
## to be 20% of xrange. For numerical reasons,
## irrespective of the value of max.window, we don't need to extend
## the smoothing window to more than 4*bandwidth.  Type is a
## categorical variable that defines the type of kernel smoothing: its
## default value "Euclidean" is the commonly used kernel smoother; the
## other value, "circular", glues the two bounds of the domain
## together therefore is more efficient for periodic (circular)
## functions.

## this is the *engine* of the fixed kernel smoother
ks.fixed <- function(x, y, bandwidth, smooth.window, type=c("Euclidean", "circular")) {
  ## match.arg gives warning when the input is not one of the
  ## declared values
  type <- match.arg(type)
  N <- length(x); yhat <- rep(0, N)
  if (type=="Euclidean") {
    for (i in 1:N){
      ## y.within.window.R contains all values to be smoothed to the
      ## right of y[j] ; dist.R are the distances (in x) from these
      ## neighboring ys to y[i]
      y.within.window.R <- c(); dist.R <- c()
      j <- 0; dj <- 0
      while (dj<=smooth.window & j<=N-i) {
        dj <- abs(x[i+j]-x[i]); dist.R <- c(dist.R, dj)
        y.within.window.R <- c(y.within.window.R, y[i+j])
        j <- j+1
      }
      ## remove the last entry if dj>smooth.window
      if (dj>smooth.window) {
        y.within.window.R <- y.within.window.R[-j]
        dist.R <- dist.R[-j]
      }
      ## now the left side
      y.within.window.L <- c(); dist.L <- c()
      j <- -1; dj <- 0
      while (dj<=smooth.window & j>=1-i) {
        dj <- abs(x[i+j]-x[i]); dist.L <- c(dist.L, dj)
        y.within.window.L <- c(y.within.window.L, y[i+j])
        j <- j-1
      }
      ## remove the last entry if dj>smooth.window
      if (dj>smooth.window) {
        y.within.window.L <- y.within.window.L[j+1]
        dist.L <- dist.L[j+1]
      }
      ## now compute the weighted average
      y.within.window <- c(y.within.window.L, y.within.window.R)
      dd <- c(dist.L, dist.R)
      w <- dnorm(dd, sd=bandwidth)
      yhat[i] <- sum(y.within.window*w)/sum(w)
    }
  } else if (type=="circular") {
    for (i in 1:N){
      ## y.within.window.R contains all values to be smoothed to the
      ## right of y[j] ; dist.R are the distances (in x) from these
      ## neighboring ys to y[i]
      y.within.window.R <- c(); dist.R <- c()
      j <- 0; dj <- 0
      while (dj<=smooth.window & j<=N-i) {
        dj <- pmin(abs(x[i+j]-x[i]), x[length(x)]- x[1] - abs(x[i+j]-x[i])) 
        dist.R <- c(dist.R, dj)
        y.within.window.R <- c(y.within.window.R, y[i+j])
        j <- j+1
      }
      ## remove the last entry if dj>smooth.window
      if (dj>smooth.window) {
        y.within.window.R <- y.within.window.R[-j]
        dist.R <- dist.R[-j]
      }
      ## now the left side
      y.within.window.L <- c(); dist.L <- c()
      j <- -1; dj <- 0
      while (dj<=smooth.window & j>=1-i) {
        dj <- pmin(abs(x[i+j]-x[i]), x[length(x)]- x[1] - abs(x[i+j]-x[i]))
        dist.L <- c(dist.L, dj)
        y.within.window.L <- c(y.within.window.L, y[i+j])
        j <- j-1
      }
      ## remove the last entry if dj>smooth.window
      if (dj>smooth.window) {
        y.within.window.L <- y.within.window.L[j+1]
        dist.L <- dist.L[j+1]
      }
      ## now compute the weighted average
      y.within.window <- c(y.within.window.L, y.within.window.R)
      dd <- c(dist.L, dist.R)
      w <- dnorm(dd, sd=bandwidth)
      yhat[i] <- sum(y.within.window*w)/sum(w)
    }
  } else {
    stop("Only two type of smoothers are currently implemented: Euclidean and circular.")
  }
  return(yhat)
}


## I decide to implement the engine of the variable bandwidth
## smoothing in a separate function, because it may be significantly
## slower than the fixed bandwidth smoother, due to the fact that many
## tests need to be performed inside of the loop. Note that both
## bandwidth and smooth.window are vectors, not a single number. They
## are assumed to have been "processed" in a way that very small
## bandwidth are truncated to zero and very large values into Inf.
ks.variable <- function(x, y, bandwidth, smooth.window, type=c("Euclidean", "circular")) {
  ## match.arg gives warning when the input is not one of the
  ## declared values
  type <- match.arg(type)
  N <- length(x); yhat <- rep(0, N)
  if (type=="Euclidean") {
    for (i in 1:N){
      if (bandwidth[i]==0){
        ## in this special case, no smoothing is needed
        yhat[i] <- y[i]
      } else { #smoothing is required
        ## y.within.window.R contains all values to be smoothed to the
        ## right of y[j] ; dist.R are the distances (in x) from these
        ## neighboring ys to y[i]
        y.within.window.R <- c(); dist.R <- c()
        j <- 0; dj <- 0
        while (dj<=smooth.window[i] & j<=N-i) {
          dj <- abs(x[i+j]-x[i]); dist.R <- c(dist.R, dj)
          y.within.window.R <- c(y.within.window.R, y[i+j])
          j <- j+1
        }
        ## remove the last entry if dj>smooth.window[i]
        if (dj>smooth.window[i]) {
          y.within.window.R <- y.within.window.R[-j]
          dist.R <- dist.R[-j]
        }
        ## now the left side
        y.within.window.L <- c(); dist.L <- c()
        j <- -1; dj <- 0
        while (dj<=smooth.window[i] & j>=1-i) {
          dj <- abs(x[i+j]-x[i]); dist.L <- c(dist.L, dj)
          y.within.window.L <- c(y.within.window.L, y[i+j])
          j <- j-1
        }
        ## remove the last entry if dj>smooth.window[i]
        if (dj>smooth.window[i]) {
          y.within.window.L <- y.within.window.L[j+1]
          dist.L <- dist.L[j+1]
        }
        y.within.window <- c(y.within.window.L, y.within.window.R)
        ## now compute the weighted average
        if (bandwidth[i]==Inf) {
          ## use simple average over the smoothing window
          yhat[i] <- mean(y.within.window)
        } else {
          ## this is the typical case
          dd <- c(dist.L, dist.R)
          w <- dnorm(dd, sd=bandwidth[i])
          yhat[i] <- sum(y.within.window*w)/sum(w)
        }
      }
    }
  } else if (type=="circular") {
    for (i in 1:N){
      if (bandwidth[i]==0){
        ## in this special case, no smoothing is needed
        yhat[i] <- y[i]
      } else { #smoothing is required
        ## y.within.window.R contains all values to be smoothed to the
        ## right of y[j] ; dist.R are the distances (in x) from these
        ## neighboring ys to y[i]
        y.within.window.R <- c(); dist.R <- c()
        j <- 0; dj <- 0
        while (dj<=smooth.window[i] & j<=N-i) {
          dj <- pmin(abs(x[i+j]-x[i]), x[length(x)]- x[1] - abs(x[i+j]-x[i])) 
          dist.R <- c(dist.R, dj)
          y.within.window.R <- c(y.within.window.R, y[i+j])
          j <- j+1
        }
        ## remove the last entry if dj>smooth.window[i]
        if (dj>smooth.window[i]) {
          y.within.window.R <- y.within.window.R[-j]
          dist.R <- dist.R[-j]
        }
        ## now the left side
        y.within.window.L <- c(); dist.L <- c()
        j <- -1; dj <- 0
        while (dj<=smooth.window[i] & j>=1-i) {
          dj <- pmin(abs(x[i+j]-x[i]), x[length(x)]- x[1] - abs(x[i+j]-x[i]))
          dist.L <- c(dist.L, dj)
          y.within.window.L <- c(y.within.window.L, y[i+j])
          j <- j-1
        }
        ## remove the last entry if dj>smooth.window[i]
        if (dj>smooth.window[i]) {
          y.within.window.L <- y.within.window.L[j+1]
          dist.L <- dist.L[j+1]
        }
        y.within.window <- c(y.within.window.L, y.within.window.R)
        ## now compute the weighted average
        if (bandwidth[i]==Inf) {
          ## use simple average over the smoothing window
          yhat[i] <- mean(y.within.window)
        } else {
          ## this is the typical case
          dd <- c(dist.L, dist.R)
          w <- dnorm(dd, sd=bandwidth[i])
          yhat[i] <- sum(y.within.window*w)/sum(w)
        }
      }
    }
  } else {
    stop("Only two type of smoothers are currently implemented: Euclidean and circular.")
  }
  return(yhat)
}




## this is the main wrapper
ks <- function(x, y, bandwidth, min.bandwidth="auto", max.bandwidth="auto", max.window="auto", type=c("Euclidean", "circular")){
  type <- match.arg(type)
  N <- length(x); xsteps <- diff(x); xrange <- max(x)-min(x)
  if (min.bandwidth=="auto") min.bandwidth <- min(xsteps)
  if (max.bandwidth=="auto") max.bandwidth <- xrange*10  # initially 10*xrange
  if (max.window=="auto") max.window <- xrange/0.1  # initially xrange/5
  ## classify the smoothing strategies based on bandwidth
  if (length(bandwidth)==1){ #fixed bandwidth
    if (bandwidth <= min.bandwidth) {
      ## we don't have to do smoothing at all
      yhat <- y
    } else if (bandwidth >= max.bandwidth) {
      ## just return the simple sample mean
      yhat <- rep(mean(y), N)
    } else { #this is the general case
      smooth.window <- min(4*bandwidth, max.window) # initially 4*bandwidth
      yhat <- ks.fixed(x, y, bandwidth, smooth.window, type=type)
    }
  } else if (length(bandwidth)!=N){
    stop("For variable bandwidth smoothing, length(bandwidth) must match that of x.")
  } else {
    ## process bandwith and smooth.window for ks.variable()
    bandwidth <- replace(bandwidth, bandwidth<=min.bandwidth, 0)
    bandwidth <- replace(bandwidth, bandwidth>=max.bandwidth, Inf) 
    smooth.window <- pmin(4*bandwidth, max.window) # initially 4*bandwidth
    yhat <- ks.variable(x, y, bandwidth, smooth.window)
  }
  return(yhat)
}
