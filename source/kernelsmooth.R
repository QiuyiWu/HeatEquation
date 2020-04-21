#--------------------- Updated Version --------------------------

kern = function(x) exp(-x^2/2) # gaussion kernel

kernsmooth = function(x,y,band){
  kij <- outer(x, x, function(x,y) kern((x-y)/(band*0.3706506)))
  s <- kij/rowSums(kij)
  return(s %*% y)    
}

               
# Example               
x <- seq(-1,1,0.1)
y <- x^2 + rnorm(length(x), 0, 0.1)                             
yhat = kernsmooth(x,y,2)      
plot(x,yhat, col=2)             
               
               

#--------------------- Old Version --------------------------

#the function gauss.kernel takes two value: a vector y and an index i

gauss.kernel = function(x, y,ind,b){
  k = exp(-((x[ind] - x)^2)/(2*b))
  yhat = sum(k*y)/sum(k)
  return(yhat)
}

gaussian.smooth = function(x,y,b){
  yhat = vector(length = length(y))
  for(i in 1:length(y)){
    yhat[i] = gauss.kernel(x,y,i,b)
  }
  return(yhat)
}
#--------------------------------------------------------------- 
             
             
#example
set.seed(1)
x <- seq(-1,1,0.1)
y <- x^2 + rnorm(length(x), 0, 0.1)
plot(x,y)
plot(x, gaussian.smooth(x,y,0.01))
plot(x, gaussian.smooth(x,y,0.2))
