

circular.dis = function(x1, xn, x, const){
  start_point = x1 
  end_point = xn
  if(x <= const){
    d1 = const - x
    d2 = (x - start_point) + (end_point - const)
    
    return(min(d1^2, d2^2))
  }else
    d1 = const - x
    d2 = (const - start_point) + (end_point - x)
    return(min(d1^2, d2^2))
}

#testing the circular distance

x = seq(1,11,1)
circular.dis(1, 10, 1,10)

#apply circular.dis for a vector
circular.dis.vec = function(x, const){
  y = rep(NA, length(x))
  for (i in 1:length(x)) {
    y[i] = circular.dis(x[1], x[length(x)],x[i],const)
   
  }
  return(y)
}

#test circular.dis.vect

z = circular.dis.vec(x,7)

gauss.kernel = function(x, y,ind,b){
 # k = exp(-((x[ind] - x)^2)/(2*b))*(1/sqrt(2*pi*b))
  k = exp(-circular.dis.vec(x,x[ind])/(2*b)) 
  #print(k)
  yhat = sum(k*y)/sum(k)
  return(yhat)
  
}



gaussian.smooth = function(x,y,b){
  yhat = vector(length = length(y))
  for(i in 1:length(y)){
    yhat[i] = gauss.kernel(x,y,i,b[i])
  }
  return(yhat)
}

################################################################

set.seed(1)
x = seq(-2*pi + 1, 2*pi + 1, pi/32)
u =  sin(x)

sd = 1
var = sd^2
y = u + rnorm(length(u),0,sd)
n = length(y)

u2 = -sin(x)# #second derivative of u
b = (var/(2*sqrt(pi)*u2^2))^(2/5) #calculate optimal bandwidth

plot(x, b, ylim = c(-0,10))
plot(x, b)
plot(x, (u2)^2)



plot(x,y, main= "Standard Deviation 100")
lines(x, u, col = "blue")
lines(x, gaussian.smooth(x, y, b), col = "red")
lines(x,gaussian.smooth(x,y, rep(0.1, length(x))), col = "green")
legend("top", legend = c("Underlying curve", "Optimal Curve"), col = c("blue", "red"), lty = 1)

SSE = (gaussian.smooth(x,y, b) - u)^2
MSE = mean(SSE)

#MSE5 = mean(SSE[c(1:390)])


bandw = 10^(seq(-8,8,1))
MSE_b = NULL
#bandw
#(gaussian.smooth(x,y,bandw[1]) - u)^2

#MSE_b5 = NULL
for (i in 1:length(bandw)) {
  SSE_vec  = (gaussian.smooth(x,y, rep(bandw[i], length(x)) ) - u)^2
  # print(SSE_vec)
  MSE_b[i] = mean(SSE_vec)
  # print(MSE_b[i])
  # MSE_b5[i] = mean(SSE_vec[1:390])
  
}
#####################################################################
#The optimal bandwidth does not have the lowest MSE. This is due to the tail behavior????

plot(seq(-8,8,1), MSE_b, ylab = "MSE", xlab = "bandwidth (10^x)", ylim = c(0, max(MSE_b)), xlim = c(-8,8), main = "MSE for different bandwidth value" )
abline(h = MSE, col = "red") #MSE of optimal bandwidth
legend("topright", legend = c("Optimal bandwidth"), col = c("red"), lty = 1)

#Since bandwidth value 10^-1 has the lowest MSE, we are going to plot and compare this with the optimal curve. 

plot(x,y, main= "Standard Deviation 100")
lines(x, u, col = "blue")
lines(x, gaussian.smooth(x, y, b), col = "red")
lines(x,gaussian.smooth(x,y, rep(0.1, length(x))), col = "green")
legend("top", legend = c("Underlying curve", "Optimal Curve", "Bandwidth of 0.1"), col = c("blue", "red", "green"), lty = 1)

min(MSE_b)
MSE
