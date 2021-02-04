# Variance with sigma


#samplesize <-c(251,501,751,1001)  # initially 501

## Note that you shouldn't use "sd" as an object; it *replaces* an
## important system function sd()!!
K = c(0.1,5,8,15)
theory.N = variance.N = array(dim = 4)
for (j in 1:4) {
  #N = samplesize[j]
  N=1001
  sd1 <- 0.5; var1 <- sd1^2  # sd initially 1
  ## an evenly-spaced grid of x defined 
  x <- seq(-pi, pi, length.out=N); dx <- range(x)/(N-1)
  L = diff(range(x))
  ## the true function to be estimated
  u <- sin(x)
  ## its second order derivative
  D2u <- -sin(x)
  ## the optimal bandwidth
  #b = (var1/(2*sqrt(pi)*D2u^2))^(2/5)
  b = 0.2
  
  simk2 <- function(y,x,u,k,b){
    yhat = ks(x, y, k,bandwidth=b, type = "circular")
    sse = (yhat - u)^2
    return(sse)
  }
  
  set.seed(123)
  y = array(dim = c(500, length(x)))
  yhat = array(dim = c(500, length(x)))
  for (m in 1:500) {
    y[m,] <- u + sd1*rnorm(N)
    #SE2[m,] = simk2(y,x,u,5,b)
    yhat[m,] = ks(x, y[m,], k=K[j], bandwidth=b, type = "circular")
    print(m)
  }
  
  variance = colMeans((yhat-matrix(rep(colMeans(yhat),500), nrow = 500, byrow = T ))^2)

  theoryMSE = NULL
  for (i in 1:N) {
    #theoryMSE[i] = sd1^2*L/((N-1)*sqrt(pi*b))
    theoryMSE[i] = sd1^2*L/((N-1)*2*sqrt(pi)*b)
  }
  
  theory.N[j] = mean(theoryMSE)
  variance.N[j] <- mean(variance)
  #variance.N[j] <- variance[30]
}



plot(K, theory.N, pch=20, ylim = c(0, max(theory.N, variance.N)+0.001), main="Variance", col = 1)
points(K, variance.N, col=2, pch=20 )
legend("topright", c("theory", "emperical"), pch=20,col=1:2)



temp = cbind(samplesize, theory.N, variance.N)

