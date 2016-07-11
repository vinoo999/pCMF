## test E[log(Z!)] with Z ~ Binom(n,p)

lgamZ = function(n,p) {
    return(sum(sapply(0:n, function(ind) return(log(factorial(ind)) * choose(n, ind) * p^ind * (1-p)^(n-ind)))))
}

lgamZ(10, 0.1)

lgamZ(37, 0.11)

lgamZ(100, 0.1)

lgamZ(170, 0.9)



pp = 0.8

res = t(sapply(1:150, function(n) {
    return(c(lgamZ(n, pp), lgamma(n*pp+1), log(factorial(round(n*pp)))))
}))

res
