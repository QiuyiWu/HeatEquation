set.seed(2569)

### Sim1: Estimating a Quadratic Function 

## specify parameters that generates data in Sim1
N <- 501
## Note that you shouldn't use "sd" as an object; it *replaces* an
## important system function sd()!!
sd1 <- 50; var1 <- sd1^2

## an evenly-spaced grid of x defined on [-5,5]
x <- seq(-5, 5, length.out=N); dx <- 10/(N-1)
## the true function to be estimated
u <- x^4+600
## its second order derivative
D2u <- 12*x^2
## the optimal bandwidth
b = (var1/(2*sqrt(pi)*D2u^2))^(2/5)

## the fixed bandwidths
min.logb <- -1; max.logb <- -0.7; nb <- 11
logbs <- seq(min.logb, max.logb, length.out=nb)

## number of repetitions
reps <- 50
plot(x, b, type="l", ylim=c(0,5), xlab="x", ylab="Adaptive bandwidth")

# Relationship between the smoothing bandwidth and mean MSE
MSEtab <- matrix(0, nrow=reps, ncol=nb+1)
colnames(MSEtab) <- c("Adaptive", paste0("logb=",logbs))

set.seed(123)
for (j in 1:reps){
  print(paste0("Repetition: ", j))
  y <- u + sd1*rnorm(N)
  yhat.adapt <- ks(x, y, bandwidth=b, type = "Euclidean")
  yhat.fixed <- sapply(logbs, function(lb) ks(x,y,bandwidth=10^lb, type = "Euclidean"))
  Res <- sweep(cbind(yhat.adapt,yhat.fixed), 1, u)
  MSEtab[j,] <- colMeans(Res^2)
}

MSE.mean <- colMeans(MSEtab); MSE.se <- apply(MSEtab, 2, sd)/sqrt(reps)

plot(logbs, MSE.mean[-1], type="b", ylim=range(MSE.mean),
     xlab="log-bandwidth", ylab="Mean MSE")
## +/- 2*STDERR
lines(logbs,MSE.mean[-1]+2*MSE.se[-1], lty=3)
lines(logbs,MSE.mean[-1]-2*MSE.se[-1], lty=3)
abline(h=MSE.mean[1], lwd=2, lty=2)


#par(mfrow=c(1,2))
plot(x, y, pch="."); lines(x,u)
lines(x, yhat.adapt, col=2)
lines(x, yhat.fixed[, 1], col=3)
lines(x, yhat.fixed[, 11], col=4)

plot(x, y, pch=1, xlim=c(4.2, 5), ylim=c(900, 1200)); lines(x,u)
lines(x, yhat.adapt, col=2)
lines(x, yhat.fixed[, 1], col=3)
lines(x, yhat.fixed[, 11], col=4)


############################################
### Sim2: Estimating a Periodic Function ###
############################################

## specify parameters that generates data in Sim2
N <- 501   # initially 501
## Note that you shouldn't use "sd" as an object; it *replaces* an
## important system function sd()!!
sd1 <- 0.1; var1 <- sd1^2  # sd initially 1

## an evenly-spaced grid of x defined 
x <- seq(-pi, pi, length.out=N); dx <- range(x)/(N-1)
## the true function to be estimated
u <- sin(x)+cos(x)
## its second order derivative
D2u <- -sin(x)-cos(x)
## the optimal bandwidth
b = (var1/(2*sqrt(pi)*D2u^2))^(2/5)


plot(x, b, type="l", ylim=c(0,5), xlab="x", ylab="Adaptive bandwidth")

## the fixed bandwidths
min.logb <- -3; max.logb <- 3; nb <- 11
logbs <- seq(min.logb, max.logb, length.out=nb)

## number of repetitions
reps <- 50

# Relationship between the smoothing bandwidth and mean MSE
MSEtab <- matrix(0, nrow=reps, ncol=nb+1)
colnames(MSEtab) <- c("Adaptive", paste0("logb=",logbs))
set.seed(123)
for (j in 1:reps){
  print(paste0("Repetition: ", j))
  y <- u + sd1*rnorm(N)
  yhat.adapt <- ks(x, y, bandwidth=b, type = "circular")
  yhat.fixed <- sapply(logbs, function(lb) ks(x,y,bandwidth=10^lb, type = "circular"))
  Res <- sweep(cbind(yhat.adapt,yhat.fixed), 1, u)
  MSEtab[j,] <- colMeans(Res^2)

}



MSE.mean <- colMeans(MSEtab); MSE.se <- apply(MSEtab, 2, sd)/sqrt(reps)
plot(logbs, MSE.mean[-1], type="b", ylim=range(MSE.mean),
     xlab="log-bandwidth", ylab="Mean MSE")
## +/- 2*STDERR
lines(logbs,MSE.mean[-1]+2*MSE.se[-1], lty=3)
lines(logbs,MSE.mean[-1]-2*MSE.se[-1], lty=3)
abline(h=MSE.mean[1], lwd=2, lty=2)


plot(x,Res[,1], ylim = c(-1,1),type = "l", ylab = "residual",main = "Resudial")
lines(x,Res[,6],  col = 2)
lines(x,Res[,7],  col = 3)
lines(x,Res[,8],  col = 4)
legend("top", c("Adaptive","logb=-0.6 (Opt)","logb=0","logb=0.6"), col =1:4, pch = 20 )


plot(x, y, ylim = c(-1,1), pch='.', main = "True function vs Estimated functions")
lines(x,u)
lines(x, yhat.adapt, col=2)
lines(x, yhat.fixed[, 3], col=7)
lines(x, yhat.fixed[, 4], col=8)
lines(x, yhat.fixed[, 5], col=3)
lines(x, yhat.fixed[, 6], col=4)
lines(x, yhat.fixed[, 7], col=5)
lines(x, yhat.fixed[, 8], col=6)
legend("topleft", c("true function","Adaptive","logb=-0.6 (Opt)","logb=0"), col =1:4, pch = 20 )



#### save plots #########

filename = paste('sim.plot11.eps')
setEPS()
postscript(filename,height=3,width=10)
layout(matrix(1:3,1,3, byrow = T))


MSE.mean <- colMeans(MSEtab); MSE.se <- apply(MSEtab, 2, sd)/sqrt(reps)
plot(logbs, MSE.mean[-1], type="b", ylim=range(MSE.mean),
     xlab="log-bandwidth", ylab="Mean MSE", main = "max.bandwidth = xrange/5")
## +/- 2*STDERR
lines(logbs,MSE.mean[-1]+2*MSE.se[-1], lty=3)
lines(logbs,MSE.mean[-1]-2*MSE.se[-1], lty=3)
abline(h=MSE.mean[1], lwd=2, lty=2)


plot(x,Res[,1], ylim = c(-1,1),type = "l", ylab = "residual",
     main = expression(paste("Resudial Plot, ",sigma," = 0.1, range = [",-pi,",",pi,"]")), col = 2)
lines(x,Res[,6],  col = 3)
lines(x,Res[,7],  col = 4)
legend("top", c("Adaptive","logb=-0.6 (Opt)", "logb=0"), col =2:4, pch = 20 )


plot(x, y, ylim = c(-1,1), pch='.', main = "True function vs Estimated functions")
lines(x,u)
lines(x, yhat.adapt, col=2)
lines(x, yhat.fixed[, 5], col=3)
lines(x, yhat.fixed[, 6], col=4)
legend("topleft", c("true function","Adaptive","logb=-0.6 (Opt)","logb=0"), col =1:4, pch = 20 )

 
dev.off()

#############




