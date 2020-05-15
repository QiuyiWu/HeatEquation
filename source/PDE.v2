# This code allows us to set the diffusive rate to different value for different x's

Laplacian <- function(y_k, dx, g){
  delta <- NULL
  len = length(y_k)
  # boundary condition
  delta[1] =  g[1]*(y_k[2]-y_k[1])/dx^2
  delta[len] = g[2]*(y_k[len-1]- y_k[len])/dx^2
  # laplacian
  delta[2: (len-1)] <- g[2: (len-1)]*diff(diff(y_k)/dx) /dx
  return(delta)
}


HeatEqn <- function(y0, dx,g, dt, Nt, history = TRUE){
  yhat = list()
  yhat[[1]] = y0
  for (i in 1:(Nt+1)) {
    yhat[[i+1]] = yhat[[i]] + Laplacian(yhat[[i]], dx,g) * dt
  }
  ifelse(history, return(yhat), return(yhat[[Nt+2]]))
}

# Example 
x <- seq(-1,1,length.out = 100)
y <- x^2 + 4 + rnorm(length(x),0,0.02) # real value  
hq1 <- HeatEqn(y0 = y, dx = x[1]-x[2], g = x^2, dt = 0.00001, Nt = 2000, history = F)
hq2 =  HeatEqn(y0 = y, dx = x[1]-x[2], g = rep(1, length(x)), dt = 0.00001, Nt = 2000, history = F)
#ks = kernsmooth(x, y, sqrt(2*0.00001*2000))

plot(x, y)
lines(x , hq1, col = "red")
lines(x, hq2, col = "blue")
legend("top", legend = c("g = 1", "g= x^2"), col = c("blue", "red"), lty = 1)
#lines(x,ks, col = "blue")
