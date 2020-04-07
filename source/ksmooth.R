# Example 
set.seed(1)
x <- seq(-1,1,0.01)
y <- x^2 + rnorm(length(x), 0, 1) # real value  
hq <- HeatEqn(y0 = y, dx = 0.001, g = 1, dt = 0.0000001, Nt = 3000, history = T)
tail(hq)
plot(x, y)
plot(x, hq[[2002]])
plot(ksmooth(x, y, "normal",  bandwidth = 2))
