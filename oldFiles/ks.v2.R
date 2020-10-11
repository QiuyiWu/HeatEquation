gauss.kernel = function(x, y,ind,b){
  k = exp(-((x[ind] - x)^2)/(2*b))*(1/(2*pi*b))
  yhat = sum(k*y)/sum(k)
  return(yhat)
}

gaussian.smooth = function(x,y,b){
  yhat = vector(length = length(y))
  for(i in 1:length(y)){
    yhat[i] = gauss.kernel(x,y,i,b)
  }
  return(yhat)
}

#simulate the data set
set.seed(1)
x = seq(-pi,pi,0.1)
u = x^4
sd = 0.1
var = sd^2
y = u + rnorm(length(u),0,sd)


u2 = 12*x^2 #second derivative of u
b = (var/(8*sqrt(pi)*u2^2))^(2/5) #calculate optimal bandwidth



plot(x,y)
lines(x, u)
lines(x, gaussian.smooth(x,y,0.01), col="red")
lines(x, gaussian.smooth(x,y,0.2), col = "blue")
lines(x, gaussian.smooth(x, y, b), col = "green")

MSE = sum((gaussian.smooth(x,y, b) - u)^2)

bandw = seq(0.001,0.2,0.001)
MSE_b = NULL
for (i in 1:length(bandw)) {
 
  MSE_b[i] = sum((gaussian.smooth(x,y,bandw[i]) - u)^2)
  
}

plot(bandw, MSE_b, xlab = "bandwidth", ylab = "MSE")
abline(h = MSE, col = "red") #MSE of optimal bandwidth

