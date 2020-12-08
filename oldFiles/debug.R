simk <- function(k,y){
  smooth.window = 2*pi/k
  
  dj <- pmin(abs(x-x[1]), x[length(x)]- x[1] - abs(x-x[1]))
  y.within.window <- y[dj <= smooth.window]
  #dd = dj[dj <= smooth.window]
  #w <- dnorm(dd, sd=b[1])
  
  #sse = (mean(y.within.window)- u[1])^2
  #return(sse)
  return(mean(y.within.window)- u[1])
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

(T.MSE = (D2u[1]*b[1]/2)^2+sd1^2/(2*sqrt(pi*b[1])))




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










