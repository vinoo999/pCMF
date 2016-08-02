### position and dispersion of gamma

## influence of shape
curve(dgamma(x, shape=1, rate=2), from=0, to=10, col="olivedrab")
curve(dgamma(x, shape=1, rate=1), from=0, to=10, col="red", add=TRUE)
curve(dgamma(x, shape=1, rate=0.5), from=0, to=10, col="blue", add=TRUE)
curve(dgamma(x, shape=1, rate=0.1), from=0, to=10, col="grey20", add=TRUE)

curve(dgamma(x, shape=1, rate=1), from=0, to=10, col="olivedrab")
curve(dgamma(x, shape=1, rate=2), from=0, to=10, col="red", add=TRUE)
curve(dgamma(x, shape=1, rate=5), from=0, to=10, col="blue", add=TRUE)
curve(dgamma(x, shape=1, rate=100), from=0, to=10, col="grey20", add=TRUE)

curve(dgamma(x, shape=0.5, rate=1), from=0, to=10, col="olivedrab")
curve(dgamma(x, shape=0.5, rate=2), from=0, to=10, col="red", add=TRUE)
curve(dgamma(x, shape=0.5, rate=5), from=0, to=10, col="blue", add=TRUE)
curve(dgamma(x, shape=0.5, rate=50), from=0, to=10, col="grey20", add=TRUE)

curve(dgamma(x, shape=2, rate=1), from=0, to=10, col="olivedrab")
curve(dgamma(x, shape=2, rate=2), from=0, to=10, col="red", add=TRUE)
curve(dgamma(x, shape=2, rate=5), from=0, to=10, col="blue", add=TRUE)
curve(dgamma(x, shape=2, rate=50), from=0, to=10, col="grey20", add=TRUE)

curve(dgamma(x, shape=3, rate=1), from=0, to=10, col="olivedrab")
curve(dgamma(x, shape=3, rate=2), from=0, to=10, col="red", add=TRUE)
curve(dgamma(x, shape=3, rate=5), from=0, to=10, col="blue", add=TRUE)
curve(dgamma(x, shape=3, rate=50), from=0, to=10, col="grey20", add=TRUE)

curve(dgamma(x, shape=0.5, rate=2), from=0, to=0.50, col="olivedrab")
curve(dgamma(x, shape=0.5, rate=1), from=0, to=0.50, col="red", add=TRUE)
curve(dgamma(x, shape=0.5, rate=0.5), from=0, to=0.50, col="blue", add=TRUE)
curve(dgamma(x, shape=0.5, rate=0.1), from=0, to=0.50, col="grey20", add=TRUE)

curve(dgamma(x, shape=2, rate=2), from=0, to=25, col="olivedrab")
curve(dgamma(x, shape=2, rate=1), from=0, to=25, col="red", add=TRUE)
curve(dgamma(x, shape=2, rate=0.5), from=0, to=25, col="blue", add=TRUE)
curve(dgamma(x, shape=2, rate=0.1), from=0, to=25, col="grey20", add=TRUE)

curve(dgamma(x, shape=2, rate=1), from=0, to=25, col="olivedrab")
curve(dgamma(x, shape=2, rate=2), from=0, to=25, col="red", add=TRUE)
curve(dgamma(x, shape=2, rate=5), from=0, to=25, col="blue", add=TRUE)
curve(dgamma(x, shape=2, rate=40), from=0, to=25, col="grey20", add=TRUE)

## influence of rate
curve(dgamma(x, shape=0.1, rate=1), from=0, to=10, col="olivedrab")
curve(dgamma(x, shape=0.5, rate=1), from=0, to=10, col="grey80", add=TRUE)
curve(dgamma(x, shape=1, rate=1), from=0, to=10, col="red", add=TRUE)
curve(dgamma(x, shape=2, rate=1), from=0, to=10, col="blue", add=TRUE)
curve(dgamma(x, shape=3, rate=1), from=0, to=10, col="grey20", add=TRUE)

curve(dgamma(x, shape=0.1, rate=0.5), from=0, to=10, col="olivedrab")
curve(dgamma(x, shape=0.5, rate=0.5), from=0, to=10, col="grey80", add=TRUE)
curve(dgamma(x, shape=1, rate=0.5), from=0, to=10, col="red", add=TRUE)
curve(dgamma(x, shape=2, rate=0.5), from=0, to=10, col="blue", add=TRUE)
curve(dgamma(x, shape=3, rate=0.5), from=0, to=10, col="grey20", add=TRUE)

curve(dgamma(x, shape=0.1, rate=2), from=0, to=10, col="olivedrab")
curve(dgamma(x, shape=0.5, rate=2), from=0, to=10, col="grey80", add=TRUE)
curve(dgamma(x, shape=1, rate=2), from=0, to=10, col="red", add=TRUE)
curve(dgamma(x, shape=2, rate=2), from=0, to=10, col="blue", add=TRUE)
curve(dgamma(x, shape=3, rate=2), from=0, to=10, col="grey20", add=TRUE)

curve(dgamma(x, shape=0.1, rate=10), from=0, to=10, col="olivedrab")
curve(dgamma(x, shape=0.5, rate=10), from=0, to=10, col="grey80", add=TRUE)
curve(dgamma(x, shape=1, rate=10), from=0, to=10, col="red", add=TRUE)
curve(dgamma(x, shape=2, rate=10), from=0, to=10, col="blue", add=TRUE)
curve(dgamma(x, shape=3, rate=10), from=0, to=10, col="grey20", add=TRUE)


curve(dgamma(x, shape=50, rate=1), from=0, to=100, col="olivedrab")




### link between a and b
a=2

curve(dgamma(x, shape=a/1, rate=1), from=0, to=10, col="olivedrab")
curve(dgamma(x, shape=a/2, rate=2), from=0, to=10, col="red", add=TRUE)
curve(dgamma(x, shape=a/4, rate=4), from=0, to=10, col="blue", add=TRUE)
curve(dgamma(x, shape=a/6, rate=6), from=0, to=10, col="grey20", add=TRUE)

curve(dgamma(x, shape=a/1, rate=1), from=0, to=10, col="olivedrab")
curve(dgamma(x, shape=a/0.8, rate=0.8), from=0, to=10, col="red", add=TRUE)
curve(dgamma(x, shape=a/0.6, rate=0.6), from=0, to=10, col="blue", add=TRUE)
curve(dgamma(x, shape=a/0.4, rate=0.4), from=0, to=10, col="grey20", add=TRUE)



curve(dgamma(x, shape=0.1, rate=1/0.1), from=0, to=10, col="olivedrab")
curve(dgamma(x, shape=0.5, rate=1/0.5), from=0, to=10, col="grey80", add=TRUE)
curve(dgamma(x, shape=2, rate=1/2), from=0, to=10, col="red", add=TRUE)
curve(dgamma(x, shape=4, rate=1/4), from=0, to=10, col="blue", add=TRUE)
curve(dgamma(x, shape=6, rate=1/6), from=0, to=10, col="grey20", add=TRUE)

curve(dgamma(x, shape=a/1, rate=1), from=0, to=10, col="olivedrab")
curve(dgamma(x, shape=0.5, rate=10), from=0, to=10, col="grey80", add=TRUE)
curve(dgamma(x, shape=a/0.8, rate=0.8), from=0, to=10, col="red", add=TRUE)
curve(dgamma(x, shape=a/0.6, rate=0.6), from=0, to=10, col="blue", add=TRUE)
curve(dgamma(x, shape=a/0.4, rate=0.4), from=0, to=10, col="grey20", add=TRUE)


