library(mvtnorm)
# Default setting
tr <- function(mat) sum(diag(as.matrix(mat)))
sn <-  0.05 # a constant in front of the bandwidth matrix (page 5)
A <- matrix(c(1, 2,2, 0), nrow=2) 
d <- nrow(A) 
Sigma0 <- crossprod(A)/det(crossprod(A))^(1/d) # determinant = 1 # bandwidth
Sigma <- sn*Sigma0
## X=c(0.15, 0.275)
dx <- 0.025
Ygrid <- as.matrix(expand.grid(seq(-2, 2, dx), seq(-2, 2, dx)))
n <- nrow(Ygrid) 
colnames(Ygrid) <- c("x", "y")
xx <- Ygrid[,1]; yy <- Ygrid[,2]
Omega <- diff(range(xx))*diff(range(yy)) #the space
cn <- Omega/n *(4*pi)^(-d/2)
### sigma setting (from eq: varZ)
a1 <- 1.2; a2 <- .8; a3 <- 50
sigma2 <- a3*((xx-a1)^2/2 +2*xx +(yy-a2)^2/4 +yy/2 +1) #increase the variance by a factor of a3 
Hsigma <- a3*matrix(c(1,0,0,0.5), nrow = 2)
#sigma2[idx]; Dsigma2x
## 2. the "true" ustar
b1 <- 1; b2 <- 0.7 
ustar <- (xx-b1)^2 +(yy-b2)^2/2
Hx <- matrix(c(2, 0, 0, 1), nrow=2)
u = drop(solve(Hx)%*%Dsigma2x)
v0 = drop(solve(Hx)%*%Dustarx)



#idx <- 1; X <- Ygrid[idx,] #x = 0.15 and y = 0.275
## gradient is c(x+2, y/2 +1/2)
Dsigma2x <- a3*c(Dx=(xx[idx]-a1)+2, Dy=(yy[idx]-a2)/2+1/2) 
Dustarx <- c(Dx=2*(xx[idx]-b1), Dy=(yy[idx]-b2)) 



MSE.noncentral = function(parameters, second.order = F, result = "MSE", k = 1){
  # result can return as just MSE, bias or variance; or return all of them: 
  # result = "MSE", result = "bias", result = "variance", result = "all". Default is "MSE"
  
  # k is the scale for delta parameter. Default is 1.

  sn = parameters[1]
  delta1 = parameters[2]
  delta2 = parameters[3]
  
  delta.n = k*c(delta1, delta2)
  second.term = 0
  if(second.order){
    second.term = 0.5*sum(drop(t(delta.n))%*%Hsigma%*%delta.n)
  }
  
  # MSE  
  MSE = (sum(Dustarx*delta.n) +sn/2*tr(Hx%*%Sigma0) +drop(t(delta.n)%*%Hx%*%delta.n)/2)^2 + cn*sn^(-d/2)*(sigma2[idx] + sum(Dsigma2x*delta.n) + second.term)
  # BIAS 
  bias = sum(Dustarx*delta.n) +sn/2*tr(Hx%*%Sigma0) +drop(t(delta.n)%*%Hx%*%delta.n)/2
  # VAR
  variance = cn*sn^(-d/2)*(sigma2[idx] + sum(Dsigma2x*delta.n) + second.term)
  
  ifelse(result == "all", return(list(MSE = MSE,Bias = bias, Biassq = bias^2, Var = variance)), 
         ifelse(result == "bias", return(bias), 
                ifelse(result == "variance", return(variance), return(MSE) )))
}

para = NULL
for (idx in 1:25921) {
  Dsigma2x <- a3*c(Dx=(xx[idx]-a1)+2, Dy=(yy[idx]-a2)/2+1/2) 
  Dustarx <- c(Dx=2*(xx[idx]-b1), Dy=(yy[idx]-b2)) 
  para[[idx]] = optim(par = c(0.01,-0.15,-0.1), fn = MSE.noncentral)$par
}

para.df <- data.frame(matrix(unlist(para), ncol = 3, byrow = T))
colnames(para.df) <- c("sn","delta1","delta2")
write.csv(para.df, file = "parameter df.csv", row.names = F)

# mapply(optim(par = c(0.01,-0.15,-0.1), fn = MSE.noncentral), 1:25) 



## 3D plot

epsilons <- rnorm(n, sd=sqrt(sigma2))
Y <- ustar+epsilons

Kn = NULL
for (idx in 1:n) {
  X <- Ygrid[idx,] 
  sn = para.df[idx,1]
  delta1 = para.df[idx,2]
  delta2 = para.df[idx,3]
  delta.n = c(delta1, delta2)
  Kn[idx] = dmvnorm(Ygrid[idx,], mean=X +delta.n, sigma=sn*Sigma0) 
}

uhatx = sum(Kn*Y) / sum(Kn)
mat = matrix(uhatx - ustar, nrow = 161, byrow = T)


persp(seq(-2, 2, dx), seq(-2, 2, dx),mat,zlab = "Difference",xlab="X-grid", ylab="Y-grid",
      theta = 30, phi = 15, axes=TRUE,scale=TRUE,box=TRUE, ticktype="detailed", 
      col = "orange", shade = 0.1)











##########################################################


# input is the observed values
# sn is a vector
# delta is a matrix

function(input, sn, delta){
  ustar = f1_input # using kernel smoothing to get it
  Dustarx = Df1_input # Gradient
  Hx = D2f1_input # Hessian
  
  # Noise
  sigma2 = f2_input
  Dsigma2 = Df2_input
  Hsigma2 = D2f2_input
  
  # Constant related to space
  Omega <- input_space #the space
  n = input_n
  d = input_dim
  cn <- Omega/n *(4*pi)^(-d/2)
  
  # Bandwidth
  A <- f_input # ?
  Sigma0 <- crossprod(A)/det(crossprod(A))^(1/d) # bandwidth
  Sigma <- sn*Sigma0
  
  
  # MSE
  MSE_FUN = function(parameters){
    
    sn = parameters[1]
    delta1 = parameters[2]
    delta2 = parameters[3]
    
    MSE = (sum(Dustarx*delta.n) +sn/2*tr(Hx%*%Sigma0) +drop(t(delta.n)%*%Hx%*%delta.n)/2)^2 + 
      cn*sn^(-d/2)*(sigma2[idx] + sum(Dsigma2x*delta.n) + 0.5*sum(drop(t(delta.n))%*%Hsigma%*%delta.n))
    return(MSE)
    }
  
  optim(c(0,0,0), MSE_FUN)  
  
  
  # Expected value of Y
  Kn = dmvnorm(Ygrid, mean=X+delta.n, sigma=sn*Sigma0) 
  uhatx = sum(Kn*Y) / sum(Kn)
  difference = abs(uhatx -  ustar)
  
  ifelse(difference > threshold, ustar = uhatx && run_again, stop) 
  
}




  






  
  
  
  
  
  
  

