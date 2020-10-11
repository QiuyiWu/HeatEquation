# This code deals with the issue of bandwidth -> infinity. 
# We set limit finite [Window] of the kernel smoother

gauss.kernel = function(x, y,ind,b,win,proportion){
  b = max(b, 1e-6)
  ran = (max(range(x)) - min(range(x)))*proportion

  xfix = x[ind]
  
  ub =  min(xfix + win*b, xfix + ran/2) #upperbound
  lb =  max(xfix - win*b, xfix - ran/2) #lowerbound
  
  xnew = x[x >= lb & x <= ub]
  ynew = y[x >= lb & x <= ub]
  k = exp(-((xfix - xnew)^2)/(2*b))
  yhat = sum(k*ynew)/sum(k)
  return(yhat)
}

gaussian.smooth = function(x,y,b,win,proportion){
  yhat = vector(length = length(y))
  for(i in 1:length(y)){
    yhat[i] = gauss.kernel(x,y,i,b[i],win,proportion)  # adaptive bandwidth
  }
  return(yhat)
}




set.seed(1)
win = 1
proportion = 1
x =  seq(-5,5,0.1) 
u =  x^4+600
sd = 100
var = sd^2
y = u + rnorm(length(u),0,sd)
n = length(y)
u2 = 12*x^2  #second derivative of u
b = (var/(2*sqrt(pi)*u2^2))^(2/5) #calculate optimal bandwidth


plot(x,y, main= "Standard Deviation 100",pch = 20)
lines(x, u, col = "blue")
lines(x, gaussian.smooth(x, y, b, win, proportion), col = "red")  # remvoe the first and last 5 points
legend("top", legend = c("Underlying curve", "Optimal Curve"), col = c("blue", "red"), bty = "n", lty = 1)

SSE = (gaussian.smooth(x,y,b)- u)^2     # remvoe the first and last 5 points
MSE = mean(SSE, na.rm = T)  # remvoe the NA point (when x=0)
bandw = 10^(seq(-8,8,0.5))
MSE_b = NULL
for (i in 1:length(bandw)) {
  SSE_vec  = (gaussian.smooth(x,y,rep(bandw[i], length(x))) - u)^2
  MSE_b[i] = mean(SSE_vec)
}

plot(seq(-8,8,0.5), MSE_b,pch = 20, ylab = "MSE", xlab = "bandwidth (10^x)", ylim = c(0, max(MSE_b)), xlim = c(-8,8), main = "MSE for different bandwidth value" )
abline(h = MSE, col = "red") #MSE of optimal bandwidth
legend("topright", legend = c("Optimal bandwidth"), col = c("red"), lty = 1,  bty = "n")
