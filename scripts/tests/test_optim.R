x = seq(-10,10,by=0.01)
f = function(x) {sin(x)}
plot(x,f(x),type="l")

x0 = numeric(20)
x0[1] = 4
for(i in 2:20) {x0[i] = cos(x0[i-1]) + x0[i-1] }
points(x0, f(x0), col="red")
x0


x = seq(-10,10,by=0.01)
f = function(x) {-(x-1)^2}
plot(x,f(x),type="l")

x0 = numeric(20)
x0[1] = -3
for(i in 2:20) {x0[i] = -2*x0[i-1] + 2 + x0[i-1] }
points(x0, f(x0), col="red")
x0
