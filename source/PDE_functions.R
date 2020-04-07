
Laplacian <- function(y_k, dx){
  delta <- NULL
  len = length(y_k)
  # boundary condition
  delta[1] = (y_k[2]-y_k[1])/dx^2
  delta[len] = (y_k[len-1]- y_k[len])/dx^2
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

# Example 
x <- seq(-1,1,0.01)
y <- x^2 + rnorm(length(x), 0, 0.05) # real value  
hq <- HeatEqn(y0 = y, dx = 0.1, g = 1, dt = 0.001, Nt = 2000, history = T)



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

 x <- seq(-1,1,0.01)
 y <- x^2 + rnorm(length(x), 0, 0.03) # real value  
 hq <- HeatEqn(y0 = y, dx = 0.1, g = 1, dt = 0.01, Nt = 200)



df <- data.frame(cbind(x,data.frame(matrix(unlist(hq[(1:6)*30]), ncol =  6))))
ggplot(df) + 
  geom_line(aes(y = X1, x = x, color = '1')) + 
  geom_line(aes(y = X2, x = x, color = '2'), linetype = "twodash") + 
  geom_line(aes(y = X3, x = x, color = '3'), linetype = "twodash") + 
  geom_line(aes(y = X4, x = x, color = '4'), linetype = "twodash") + 
  geom_line(aes(y = X5, x = x, color = '5'), linetype = "twodash") + 
  geom_line(aes(y = X6, x = x, color = '6'), linetype = "twodash") + 
  scale_color_discrete(name = "Y", labels = c("real value", "fitted value 1", "fitted value 2",
                                              "fitted value 3", "fitted value 4", "fitted value 5")) +
  ggtitle("One Dimensional Finite Difference Method")





