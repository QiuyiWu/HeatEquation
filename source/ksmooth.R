# Example 
set.seed(1)
x <- seq(-1,1,0.1)
y <- x^2 + rnorm(length(x), 0, 0.1) # real value  
hq <- HeatEqn(y0 = y, dx = 0.1, g = 1, dt = 0.001, Nt = 275, history = F)
tail(hq)
plot(x, y)
plot(x, hq)
plot(ksmooth(x, y, "normal",  bandwidth = 2))
