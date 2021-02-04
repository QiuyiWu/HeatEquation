N=1001
sd1 <- 0.5; var1 <- sd1^2  # sd initially 1
## an evenly-spaced grid of x defined 
x <- seq(-pi, pi, length.out=N); dx <- range(x)/(N-1)
L = diff(range(x))
## the true function to be estimated
u <- x^4
## its second order derivative
D2u <- 12*x^2
## the optimal bandwidth
#b = (var1/(2*sqrt(pi)*D2u^2))^(2/5)
b = 0.2

simk2 <- function(y,x,u,k,b){
  yhat = ks(x, y, k,bandwidth=b, type = "circular")
  sse = (yhat - u)^2
  return(sse)
}

set.seed(123)
y = array(dim = c(100, length(x)))
yhat = array(dim = c(100, length(x)))
for (m in 1:100) {
  y[m,] <- u + sd1*rnorm(N)
  #SE2[m,] = simk2(y,x,u,5,b)
  yhat[m,] = ks(x, y[m,], k=5, bandwidth=b, type = "circular")
  print(m)
}

biass = (colMeans(yhat)-u)^2



theoryMSE = NULL
for (i in 1:N) {
  theoryMSE[i] = (0.5*D2u[i]*b^2)^2
}

plot(x, theoryMSE, main="Bias (N=1001, K=5), y=x^4",col = 1,pch=19, xlim = c(-2,2),ylim = c(0,1))
points(x, biass, col=2, pch=20 )
legend("top", c("theory", "emperical"), pch=c(19,20),col=1:2, bty="n")

