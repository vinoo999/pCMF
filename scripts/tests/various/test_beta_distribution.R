### position and dispersion of gamma

## influence of shape
curve(dbeta(x, shape=1, shape2=5), from=0, to=1, col="olivedrab")
curve(dbeta(x, shape=1, shape2=3), from=0, to=1, col="green", add=TRUE)
curve(dbeta(x, shape=1, shape2=2), from=0, to=1, col="orange", add=TRUE)
curve(dbeta(x, shape=1, shape2=1), from=0, to=1, col="grey80", add=TRUE)
curve(dbeta(x, shape=1, shape2=0.8), from=0, to=1, col="red", add=TRUE)
curve(dbeta(x, shape=1, shape2=0.5), from=0, to=1, col="blue", add=TRUE)
curve(dbeta(x, shape=1, shape2=0.1), from=0, to=1, col="grey20", add=TRUE)

curve(dbeta(x, shape=5, shape2=1), from=0, to=1, col="olivedrab")
curve(dbeta(x, shape=3, shape2=1), from=0, to=1, col="green", add=TRUE)
curve(dbeta(x, shape=2, shape2=1), from=0, to=1, col="orange", add=TRUE)
curve(dbeta(x, shape=1, shape2=1), from=0, to=1, col="grey80", add=TRUE)
curve(dbeta(x, shape=0.8, shape2=1), from=0, to=1, col="red", add=TRUE)
curve(dbeta(x, shape=0.5, shape2=1), from=0, to=1, col="blue", add=TRUE)
curve(dbeta(x, shape=0.1, shape2=1), from=0, to=1, col="grey20", add=TRUE)

curve(dbeta(x, shape=0.5, shape2=2), from=0, to=0.50, col="olivedrab")
curve(dbeta(x, shape=0.5, shape2=1), from=0, to=0.50, col="red", add=TRUE)
curve(dbeta(x, shape=0.5, shape2=0.5), from=0, to=0.50, col="blue", add=TRUE)
curve(dbeta(x, shape=0.5, shape2=0.1), from=0, to=0.50, col="grey20", add=TRUE)

curve(dbeta(x, shape=2, shape2=2), from=0, to=25, col="olivedrab")
curve(dbeta(x, shape=2, shape2=1), from=0, to=25, col="red", add=TRUE)
curve(dbeta(x, shape=2, shape2=0.5), from=0, to=25, col="blue", add=TRUE)
curve(dbeta(x, shape=2, shape2=0.1), from=0, to=25, col="grey20", add=TRUE)

curve(dbeta(x, shape=2, shape2=1), from=0, to=25, col="olivedrab")
curve(dbeta(x, shape=2, shape2=2), from=0, to=25, col="red", add=TRUE)
curve(dbeta(x, shape=2, shape2=5), from=0, to=25, col="blue", add=TRUE)
curve(dbeta(x, shape=2, shape2=40), from=0, to=25, col="grey20", add=TRUE)

## influence of shape2
curve(dbeta(x, shape=0.1, shape2=1), from=0, to=10, col="olivedrab")
curve(dbeta(x, shape=0.5, shape2=1), from=0, to=10, col="grey80", add=TRUE)
curve(dbeta(x, shape=1, shape2=1), from=0, to=10, col="red", add=TRUE)
curve(dbeta(x, shape=2, shape2=1), from=0, to=10, col="blue", add=TRUE)
curve(dbeta(x, shape=3, shape2=1), from=0, to=10, col="grey20", add=TRUE)

curve(dbeta(x, shape=0.1, shape2=0.5), from=0, to=10, col="olivedrab")
curve(dbeta(x, shape=0.5, shape2=0.5), from=0, to=10, col="grey80", add=TRUE)
curve(dbeta(x, shape=1, shape2=0.5), from=0, to=10, col="red", add=TRUE)
curve(dbeta(x, shape=2, shape2=0.5), from=0, to=10, col="blue", add=TRUE)
curve(dbeta(x, shape=3, shape2=0.5), from=0, to=10, col="grey20", add=TRUE)

curve(dbeta(x, shape=0.1, shape2=2), from=0, to=10, col="olivedrab")
curve(dbeta(x, shape=0.5, shape2=2), from=0, to=10, col="grey80", add=TRUE)
curve(dbeta(x, shape=1, shape2=2), from=0, to=10, col="red", add=TRUE)
curve(dbeta(x, shape=2, shape2=2), from=0, to=10, col="blue", add=TRUE)
curve(dbeta(x, shape=3, shape2=2), from=0, to=10, col="grey20", add=TRUE)

curve(dbeta(x, shape=0.1, shape2=10), from=0, to=10, col="olivedrab")
curve(dbeta(x, shape=0.5, shape2=10), from=0, to=10, col="grey80", add=TRUE)
curve(dbeta(x, shape=1, shape2=10), from=0, to=10, col="red", add=TRUE)
curve(dbeta(x, shape=2, shape2=10), from=0, to=10, col="blue", add=TRUE)
curve(dbeta(x, shape=3, shape2=10), from=0, to=10, col="grey20", add=TRUE)


curve(dbeta(x, shape=50, shape2=1), from=0, to=100, col="olivedrab")




### link between a and b
a=2

curve(dbeta(x, shape=a/1, shape2=1), from=0, to=10, col="olivedrab")
curve(dbeta(x, shape=a/2, shape2=2), from=0, to=10, col="red", add=TRUE)
curve(dbeta(x, shape=a/4, shape2=4), from=0, to=10, col="blue", add=TRUE)
curve(dbeta(x, shape=a/6, shape2=6), from=0, to=10, col="grey20", add=TRUE)

curve(dbeta(x, shape=a/1, shape2=1), from=0, to=10, col="olivedrab")
curve(dbeta(x, shape=a/0.8, shape2=0.8), from=0, to=10, col="red", add=TRUE)
curve(dbeta(x, shape=a/0.6, shape2=0.6), from=0, to=10, col="blue", add=TRUE)
curve(dbeta(x, shape=a/0.4, shape2=0.4), from=0, to=10, col="grey20", add=TRUE)



curve(dbeta(x, shape=0.1, shape2=1/0.1), from=0, to=10, col="olivedrab")
curve(dbeta(x, shape=0.5, shape2=1/0.5), from=0, to=10, col="grey80", add=TRUE)
curve(dbeta(x, shape=2, shape2=1/2), from=0, to=10, col="red", add=TRUE)
curve(dbeta(x, shape=4, shape2=1/4), from=0, to=10, col="blue", add=TRUE)
curve(dbeta(x, shape=6, shape2=1/6), from=0, to=10, col="grey20", add=TRUE)

curve(dbeta(x, shape=a/1, shape2=1), from=0, to=10, col="olivedrab")
curve(dbeta(x, shape=0.5, shape2=10), from=0, to=10, col="grey80", add=TRUE)
curve(dbeta(x, shape=a/0.8, shape2=0.8), from=0, to=10, col="red", add=TRUE)
curve(dbeta(x, shape=a/0.6, shape2=0.6), from=0, to=10, col="blue", add=TRUE)
curve(dbeta(x, shape=a/0.4, shape2=0.4), from=0, to=10, col="grey20", add=TRUE)


