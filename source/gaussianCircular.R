

Laplacian <- function(y_k, dx){
  delta <- NULL
  len = length(y_k)
  # boundary condition
  delta[1] =(y_k[2]-y_k[1])/dx^2 
  delta[len] = (y_k[len-1]- y_k[len])/dx^2
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


kern = function(x) exp(-x^2/2) # gaussion kernel

kernsmooth = function(x,y,band){
  kij <- outer(x, x, function(x,y) kern((x-y)/(band*0.3706506))) 
  s <- kij/rowSums(kij)
  return(s %*% y)    
}



#example 
x <- seq(0,2*pi,0.1)
y <- sin(x) # real value  
hq <- HeatEqn(y0 = y, dx = 0.1, g = 1, dt = 0.001, Nt = 69, history = F)
ks = ksmooth(x,y, kernel = "normal", bandwidth = 0.3706505)
ks2 = kernsmooth(x,y,1)


plot(x, y)
lines(x , hq, col = "red")
lines(ks$x, ks$y, col = "blue")
lines(x, ks2, col = "green")
legend("topright", legend = c("Finite Diff", "Ksmooth"), col= c("red", "blue"), lty = 1)


#the relationship between the ksmooth and finite difference

(((1*0.3706506)^2)/2)/0.001
