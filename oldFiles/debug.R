simk <- function(k,y){
  smooth.window = 2*pi/k
  y.within.window.R <- c(); dist.R <- c()
  j <- 0; dj <- 0
  while (dj<=smooth.window & j<=N-i) {
    dj <- pmin(abs(x[i+j]-x[i]), x[length(x)]- x[1] - abs(x[i+j]-x[i])) 
    dist.R <- c(dist.R, dj)
    y.within.window.R <- c(y.within.window.R, y[i+j])
    j <- j+1
  }
  
  if (dj>smooth.window) {
    y.within.window.R <- y.within.window.R[-j]
    dist.R <- dist.R[-j]
  }
  
  dd = dist.R
  y.within.window = y.within.window.R
  dnorm(dd, sd=b[1])
  
  sse = (mean(y.within.window)- u[1])^2
  return(sse)
  #return(mean(y.within.window)- u[1])
}


k=5
SE = NULL
for (m in 1:1000) {
  y <- u + sd1*rnorm(N)
  SE[m] = simk(k,y)
}
(MSE = mean(SE))
(VAR = var(SE))
(D2u[1]*b[1]/2)^2+sd1^2/(2*sqrt(pi*b[1]))

(T.MSE = (D2u[126]*b/2)^2+sd1^2/(2*sqrt(pi*b)))




set.seed(123)
y <- u + 1*rnorm(N)


filename = paste('sim.plot4.eps')
setEPS()
postscript(filename,height=6,width=16)
layout(matrix(1:2,1,2, byrow = T))

plot(seq(0.1,20, 0.05),unlist(lapply(seq(0.1,20, 0.05), simk, y=y)),pch =20, 
     xlab = "k", ylab = "difference between Yhat and u", main = "sd=1")
abline(v = 2, col=2, lwd = 3)
plot(seq(-3,8, 1),unlist(lapply(10^(seq(-3,8, 1)), simk, y=y)),pch = 20,
     xlab = "10^k", ylab = "difference between Yhat and u", main = "sd=1, log k scale")
abline(v = log(2), col=2,lwd = 3)

dev.off()




set.seed(123)
y <- u + 0.1*rnorm(N)

filename = paste('sim.plot5.eps')
setEPS()
postscript(filename,height=6,width=16)
layout(matrix(1:2,1,2, byrow = T))

plot(seq(0.1,20, 0.05),unlist(lapply(seq(0.1,20, 0.05), simk, y=y)),pch = 20, 
     xlab = "k", ylab = "difference between Yhat and u", main = "sd=0.1")
abline(v = 2, col=2, lwd = 3)
plot(seq(-3,8, 1),unlist(lapply(10^(seq(-3,8, 1)), simk, y=y)),pch = 20,
     xlab = "10^k", ylab = "difference between Yhat and u", main = "sd=0.1, log k scale")
abline(v = log(2), col=2,lwd = 3)

dev.off()


(bmin1 <- b[which(b==min(b))[1]]) #126
(bmin2 <- b[which(b==min(b))[2]]) #376

b = bmin1










