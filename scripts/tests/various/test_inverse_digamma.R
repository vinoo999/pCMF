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

    return(x)
}

y = digamma(0.0000002)
y
psiInv(y)



y = c(digamma(0.478), digamma(1), digamma(1.54895158), digamma(20))
y
sapply(y, psiInv)
