
Laplacian <- function(y_k, dx){
  diff(y_k)/dx  #dx = 0.1
}

HeatEqu <- function(y0, dx, g, dt, Nt){
  yhat = array(NA, Nt)
  yhat[1] = y0
  for (i in 1:(Nt-1)) {
    yhat[i+1] = yhat[i] + Laplacian(y, dx = 0.1)[i] * 0.01
  }
  return(yhat)
}

#---- Example ------

# x <- seq(-1,1,0.01)
# y <- x^2 + rnorm(length(x), 0, 0.03) # real value  
# HeatEqu(y0 = y[1], dx = 0.1, g = 1, dt = 0.01, Nt = 200)

