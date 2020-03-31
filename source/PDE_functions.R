
Laplacian <- function(y_k, dx){
  diff(diff(y_k)/dx) /dx
}

#Laplacian = function(y, dx){
#  len = length(y)
#  delta2 = NULL
#  delta2[1] = (y[2]-y[1])/dx^2
#  delta2[len] = (y[len - 1] - y[len])/dx^2

#  for(i in 2:(len - 1)){
#    delta2[i] = (y[i+1] + y[i-1] - 2*y[i])/dx^2
# } 
#  return(delta2)
#}

HeatEqn <- function(y0, dx, g, dt, Nt){
  yhat = array(NA, Nt)
  yhat[1] = y0
  for (i in 1:(Nt-1)) {
    yhat[i+1] = yhat[i] + g* Laplacian(y, dx = 0.1)[i] * dt
  }
  return(yhat)
}


#HeatEqn = function(y0, dx, g, dt, Nt, HISTORY = TRUE){
  
#  y = list()
#  y[[1]] = y0
  
#  for (i in 2:(Nt+1)){ 
#    y[[i]] = y[[i-1]] + dt*g*Laplacian(y[[i-1]],dx)
#  }
  
#  if(HISTORY){
#    return(y)
#  }else{
#    return(y[[Nt+1]])
#  }
#}

# ------ Example -------

# x <- seq(-1,1,0.01)
# y <- x^2 + rnorm(length(x), 0, 0.03) # real value  
# HeatEqn(y0 = y[1], dx = 0.1, g = 1, dt = 0.01, Nt = 200)




