#This code deals with the issue of boundary point. 
#Note: this method works well only with seasonal data (aka. the data with sine or cosine curve)

gauss.kernel = function(x, y,ind,b){
  dv = pmin(abs(x[ind] - x), x[length(x)]- x[1] - abs(x[ind] - x))  # circular boundary (pick up the minimum distance)
  k = exp(-(dv^2)/(2*b))*(1/sqrt(2*pi*b))
  yhat = sum(k*y)/sum(k)
  return(yhat)
}

gaussian.smooth = function(x,y,b){
  yhat = vector(length = length(y))
  for(i in 1:length(y)){
    yhat[i] = gauss.kernel(x,y,i,b[i])  # adaptive bandwidth
  }
  return(yhat)
}

################################################################
####################     y = x^4    ############################
################################################################

set.seed(1)
x =  seq(-5,5,0.1) 
u =  x^4 

sd = 100
var = sd^2
y = u + rnorm(length(u),0,sd)
n = length(y)

u2 = 12*x^2  #second derivative of u
b = (var/(2*sqrt(pi)*u2^2))^(2/5) #calculate optimal bandwidth

plot(x,y, main= "Standard Deviation 100",pch = 20)
lines(x, u, col = "blue")
lines(x[-c(1:5,97:101)], gaussian.smooth(x, y, b)[-c(1:5,97:101)], col = "red")  # remvoe the first and last 5 points
legend("top", legend = c("Underlying curve", "Optimal Curve"), col = c("blue", "red"), bty = "n", lty = 1)

SSE = (gaussian.smooth(x,y, b)[-c(1:5,97:101)] - u[-c(1:5,97:101)])^2     # remvoe the first and last 5 points
MSE = mean(SSE, na.rm = T)  # remvoe the NA point (when x=0)

bandw = 10^(seq(-8,8,length.out = 101))
MSE_b = NULL
for (i in 1:length(bandw)) {
  SSE_vec  = (gaussian.smooth(x,y,rep(bandw[i], length(x)) ) - u)^2
  MSE_b[i] = mean(SSE_vec)
}

plot(seq(-8,8,length.out = 101), MSE_b,pch = 20, ylab = "MSE", xlab = "bandwidth (10^x)", ylim = c(0, max(MSE_b)), xlim = c(-8,8), main = "MSE for different bandwidth value" )
abline(h = MSE, col = "red") #MSE of optimal bandwidth
legend("topleft", legend = c("Optimal bandwidth"), col = c("red"), lty = 1,  bty = "n")



################################################################
####################   y = sin(x)   ############################
################################################################


set.seed(1)
x = seq(-2*pi, 2*pi, pi/32)
u =  sin(x)

sd = 1
var = sd^2
y = u + rnorm(length(u),0,sd)
n = length(y)

u2 = -sin(x)   #second derivative of u
b = (var/(2*sqrt(pi)*u2^2))^(2/5) #calculate optimal bandwidth


plot(x,y, main= "Standard Deviation 100", pch = 20)
lines(x, u, col = "blue")
lines(x, gaussian.smooth(x, y, b), col = "red")
legend("top", legend = c("Underlying curve", "Optimal Curve"), bty = 'n', col = c("blue", "red"), lty = 1)

SSE = (gaussian.smooth(x,y, b) - u)^2
MSE = mean(SSE, na.rm = T)

bandw = 10^(seq(-8,8,length.out = 129))
MSE_b = NULL
for (i in 1:length(bandw)) {
  SSE_vec  = (gaussian.smooth(x,y,rep(bandw[i], length(x))) - u)^2
  MSE_b[i] = mean(SSE_vec)
}

plot(seq(-8,8,length.out = 129), MSE_b,pch = 20, ylab = "MSE", xlab = "bandwidth (10^x)", ylim = c(0, max(MSE_b)), xlim = c(-8,8), main = "MSE for different bandwidth value" )
abline(h = MSE, col = "red") #MSE of optimal bandwidth
legend("topright", legend = c("Optimal bandwidth"), col = c("red"), lty = 1,  bty = "n")


