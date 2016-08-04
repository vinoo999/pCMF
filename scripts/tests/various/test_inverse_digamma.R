### test inverse digamma function

psiInv = function(y, iter=5) {
    x0 = 0
    if( y >= -2.22) {
        x0 = exp(y) + 0.5
    } else {
        x0 = -1/(y-digamma(1))
    }
    x = numeric(iter+1)
    x[1] = x0
    for(i in 1:iter) {
        x[i+1] = x[i] - (digamma(x[i]) - y)/trigamma(x[i])
    }

    return(tail(x,1))
}

y = digamma(0.0000002)
y
psiInv(y)

digammaInv = function(y) {
    return(sapply(y, psiInv))
}



y = c(digamma(0.478), digamma(1), digamma(1.54895158), digamma(20))
y
sapply(y, psiInv)

curve(digamma(x), from=0, to=10)
curve(digammaInv(x), from=0, to=10)


#### test system resolution with inverse of digamma
K=4

a1 = numeric(100)
a2 = numeric(100)

a1[1] = 5
a2[1] = 4

for(i in 1:99) {
    a1[i+1] = psiInv(log(a2[i]) + 3)
    a2[i+1] = a1[i+1] / 10
}

at1 = numeric(100)
at2 = numeric(100)

at1[1] = a1[1] / sqrt(K)
at2[1] = a2[1] / sqrt(K)

for(i in 1:99) {
    at1[i+1] = psiInv(log(sqrt(K)*at2[i]) + 3)
    at2[i+1] = at1[i+1] / (sqrt(K) * 10)
}

a1[100]
at1[100] * sqrt(K)

a2[100]
at2[100] * sqrt(K)

