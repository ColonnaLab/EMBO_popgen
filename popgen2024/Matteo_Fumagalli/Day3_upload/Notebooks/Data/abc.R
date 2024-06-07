

simulateElephants <- function(param) rpois(n=1, lambda=param)

simulateWaterTemp <- function(param) rnorm(n=1, mean=param, sd=sqrt(10)) * rbeta(n=1, shape1=param^2, shape2=1/param)








