# simularion for 1000 times
# plot MSE vs bandwidth of fix/optimal bandwidth
# plot MSE vs x at every point

##################### Function ##################################
# Gaussian kernel without window 
gauss.kernel2 = function(x, y,ind,b){
  k = exp(-((x[ind] - x)^2)/(2*b))*(1/sqrt(2*pi*b))
  yhat = sum(k*y)/sum(k)
  return(yhat)
}
gaussian.smooth2 = function(x,y,b){
  yhat = vector(length = length(y))
  for(i in 1:length(y)){
    yhat[i] = gauss.kernel2(x,y,i,b)
  }
  return(yhat)
}

# Gaussian kernel with window 
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


##################### Example ################################

win = 2
proportion = 0.2
x =  seq(-5,5,length.out = 200) 
u =  x^4+600
sd = 100
var = sd^2
#y = u + rnorm(length(u),0,sd)
n = length(x)
u2 = 12*x^2  #second derivative of u
b = (var/(2*sqrt(pi)*u2^2))^(2/5) #calculate optimal bandwidth

plot(x,b, ylim = c(0,10), ylab = "bandwidth", main = "Value of bandwidth vs. x value")

MSE_b = NULL
bandw = 10^(seq(-8,8,1))

SSE.true = matrix(rep(NA, 1000*n),nrow = 1000)
MSE_b = matrix(rep(NA, 1000*length(bandw)), nrow = 1000)
MSE = NULL
for (j in 1:1000) {
  set.seed(j)
  y = u + rnorm(length(u),0,sd)
  SSE = (gaussian.smooth(x,y,b, win, proportion ) - u)^2
  MSE[j] = mean(SSE)
  SSE.true[j,] = SSE
  for (i in 1:length(bandw)) {
    SSE_vec  = (gaussian.smooth2(x,y, rep( bandw[i], length(b))) - u)^2
    MSE_b[j,i] = mean(SSE_vec)
  }
}


plot(seq(-8,8,1), colMeans(MSE_b), type = "p", ylab = "MSE", xlab = "bandwidth (10^x)", ylim = c(0, max( colMeans(MSE_b) )), xlim = c(-8,8), pch=20, main = "MSE for different bandwidth value" )
abline(h = mean(MSE), col = "red") #MSE of optimal bandwidth
legend("topleft", pch = c(95, 20), legend = c("Optimal bw = 719.38", "min(fixed bw) = 816.93"), col = c("red", "black"))

plot(x, y = colMeans(SSE.true), pch = 20, ylab = "Estimated MSE", main = "Estimation of MSE for each point")

