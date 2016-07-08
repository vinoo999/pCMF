curve((1-x)*log(x/(1-x)), from=0, to=1)
curve(expit(log(x*0.001/(1-x))), from=0.99, to=1)
curve(exp(x)/(1+exp(x)), from=-10, to=10)


logit = function(x) {
    return(log(x/(1-x)))
}

expit = function(x) {
    exp(x)/(1+exp(x))
}

p = 0.8
logit(p)

expit(0.9)



curve(logit(x), from=0, to=1)
curve(expit(x), from=-5, to=5)
