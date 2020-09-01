gauss.kernel = function(x, y,ind,b){
  dv = pmin(abs(x[ind] - x), x[length(x)]- x[1] - abs(x[ind] - x))
  k = exp(-(dv^2)/(2*b))*(1/sqrt(2*pi*b))
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
