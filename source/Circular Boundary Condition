#This code deals with the issue of boundary point. 
#Note: this method works well only with seasonal data (aka. the data with sine or cosine curve)

Laplacian <- function(y_k, dx){
  delta <- NULL
  len = length(y_k)
  # boundary condition
  delta[1] = ((y_k[2]-y_k[1])/dx - (y_k[1] - y_k[len])/dx)/dx  #the change in the first point is calculated using the last point in the array. 
  delta[len] = ((y_k[1] - y_k[len])/dx - (y_k[len]- y_k[len-1])/dx)/dx #the change in the last point is calculated using the first point.
  # laplacian
  delta[2: (len-1)] <- diff(diff(y_k)/dx) /dx
  return(delta)
}


HeatEqn <- function(y0, dx, g, dt, Nt, history = TRUE){
  yhat = list()
  yhat[[1]] = y0
  for (i in 1:(Nt+1)) {
    yhat[[i+1]] = yhat[[i]] + g* Laplacian(yhat[[i]], dx) * dt
  }
  ifelse(history, return(yhat), return(yhat[[Nt+2]]))
}

