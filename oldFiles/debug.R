setwd("~/Desktop/Research/UR/HeatEquation/latex/pic")
a=11

simk <- function(k,y,a){
  smooth.window = 2*pi/k
  
  dj <- pmin(abs(x-x[a]), x[length(x)]- x[a] - abs(x-x[a]))
  y.within.window <- y[dj <= smooth.window]
  
  
  #dd = dj[dj <= smooth.window]
  #w <- dnorm(dd, sd=b[a])
  
  sse = (mean(y.within.window)- u[a])^2
  return(sse)
  #return(mean(y.within.window)- u[a])
}

# k=5
# SE = NULL
# for (m in 1:1000) {
#   y <- u + sd1*rnorm(N)
#   SE[m] = simk(k,y)
# }
# (MSE = mean(SE))
# # (VAR = var(SE))
# (D2u[a]*b[a]/2)^2+sd1^2/(2*sqrt(pi*b[a]))
# 

k=25
SE = array(dim = c(1000, length(x)))
for (m in 1:1000) {
  for (a in 1:length(x)) {
    y <- u + sd1*rnorm(N)
    SE[m,a] = simk(k,y,a)
  }
}



theoryMSE = NULL
for (i in 1:501) {
  theoryMSE[i] = (D2u[i]*b[i]/2)^2+sd1^2/(2*sqrt(pi*b[i]))
}
diff = colMeans(SE) - theoryMSE
plot(x, diff)



####################################################################################
N <-1001   # initially 501
## Note that you shouldn't use "sd" as an object; it *replaces* an
## important system function sd()!!
sd1 <- 0.1; var1 <- sd1^2  # sd initially 1

## an evenly-spaced grid of x defined 
x <- seq(-pi, pi, length.out=N); dx <- range(x)/(N-1)
## the true function to be estimated
u <- sin(x)
## its second order derivative
D2u <- -sin(x)
## the optimal bandwidth
#b = (var1/(2*sqrt(pi)*D2u^2))^(2/5)
b = rep(0.2, N)

simk <- function(y,x,u,k,b){
  yhat = ks(x, y, k,bandwidth=b, type = "Euclidean")
  sse = (yhat - u)^2
  return(sse)
}

simk2 <- function(y,x,u,k,b){
  yhat = ks(x, y, k,bandwidth=b, type = "circular")
  sse = (yhat - u)^2
  return(sse)
}


SE = array(dim = c(50, length(x)))
for (m in 1:50) {
    y <- u + sd1*rnorm(N)
    SE[m,] = simk(y,x,u,5,b)
    print(m)
}

SE2 = array(dim = c(500, length(x)))
y = array(dim = c(500, length(x)))
yhat = array(dim = c(500, length(x)))
for (m in 1:500) {
  y[m,] <- u + sd1*rnorm(N)
  #SE2[m,] = simk2(y,x,u,5,b)
  yhat[m,] = ks(x, y[m,], bandwidth=b, type = "circular")
  print(m)
}

variance = colMeans((yhat-matrix(rep(colMeans(yhat),500), nrow = 500, byrow = T ))^2)
biass = (colMeans(yhat)-u)^2

theoryMSE = NULL
theoryMSE2 = NULL
theoryMSE3 = NULL
for (i in 1:N) {
  theoryMSE[i] = (0.5*D2u[i]*b[i]^2)^2 #+sd1^2/(2*sqrt(pi*b[i]))#+
  theoryMSE2[i] =   20*(D2u[i]*b[i]/(2*N))^2 #+sd1^2/(2*N*sqrt(pi*b[i]))# +
  theoryMSE3[i] =   20*(D2u[i]^2*b[i]/(2*sqrt(N)))^2 #+sd1^2/(2*sqrt(pi*b[i]*N))# +
}
diff = colMeans(SE) - theoryMSE
diff2 = colMeans(SE2) - theoryMSE
#apply(SE, 2, sd)

plot(x[-c(1:20,480:501)], diff[-c(1:20,480:501)], type = 'l')




filename = paste('sim.diff1.eps')
setEPS()
postscript(filename,height=6,width=15)
layout(matrix(1:2,1,2, byrow = T))

plot(x,diff, pch = 20, col=4,main = "MSE Diff (k=5, N=501, N_sim=500)", ylab = 'MSE difference')
points(x, diff2, col=2, pch = 20,title = "Euclidean")
legend('top', c("circular",'euclidean'), col = c(2,4),pch = 20, bty='n')

plot(x,diff2, pch = 20, col=2,main = "Zoom-In version", ylab = 'MSE difference')
points(x, diff, pch=20,col=4, title = "Euclidean")
legend('top', c("circular",'euclidean'), col = c(2,4),pch = 20, bty='n')
dev.off()


filename = paste('sim.mse3.eps')
setEPS()
postscript(filename,height=4,width=10)
layout(matrix(1:2,1,2, byrow = T))


plot(x, theoryMSE,pch='.', main="Variance" , ylim = c(-0.005, 0.007), ylab = "variance")
points(x, theoryMSE2,pch='.', main="Their Theoretical MSE", col= 2)
points(x, theoryMSE3,pch='.', main="Their Theoretical MSE", col= 3 )
points(x, variance, pch='.', col= 4)
legend("bottom", c("var1","var2","var3","emperical"),pch=20, col = 1:4, bty='n')

plot(x, theoryMSE2,pch='.', main="Zoom-In Variance", col= 2, ylim = c(-0.0001, 0.001), ylab = "variance")
points(x, theoryMSE3,pch='.', main="Their Theoretical MSE", col= 3 )
points(x, variance, pch='.', col= 4)
legend("top", c("var2","var3","emperical"), pch=20, col = 2:4, bty='n')


dev.off()






filename = paste('sim.mse4.eps')
setEPS()
postscript(filename,height=4,width=10)
layout(matrix(1:2,1,2, byrow = T))


plot(x, theoryMSE,pch='.', main="Bias^2" , ylim = c(-0.01, 0.01), ylab="bias^2")
points(x, theoryMSE2,pch='.', main="Their Theoretical MSE", col= 2)
points(x, theoryMSE3,pch='.', main="Their Theoretical MSE", col= 3 )
points(x, biass, pch='.', col= 4)
legend("bottom", c("bs1","bs2","bs3","emperical"),pch=20,col = 1:4, bty='n')

plot(x, theoryMSE2,pch='.', main="Zoom-In Bias^2", col= 2, ylim = c(-0.0001, 0.001), ylab="bias^2")
points(x, theoryMSE3,pch='.', main="Their Theoretical MSE", col= 3 )
points(x, biass, pch='.', col= 4)
legend("top", c("bs2","bs3","emperical"), pch=20, col = 2:4, bty='n')


dev.off()





####################################################################################













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










