curve(log(x), from=0, to=1)

curve(lgamma(x+1), from=0, to=10)
curve(digamma(x), from=0, to=10)
curve(trigamma(x), from=0, to=10)
abline(a=0,b=1)



res = t(sapply(1:100000, function(n) {
    a1 = lgamma(n)
    a2 = n*log(n) - n + 0.5*log(2*pi*n)
    return(c(a1, a2, abs(a1-a2), round(100*abs(a1-a2)/a2,2)))
}))

res
