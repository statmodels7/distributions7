# Weibull Random Number Generator

Generates random numbers by inverse transform, which is exact here
because the quantile function is elementary.

## Arguments

- distrib:

  A `Weibull1Distrib` object.

- n:

  Number of observations to generate.

- theta:

  A list containing the parameters `mu` and `sigma`.

## Value

A numeric vector of random draws.

## See also

[`weibull1_distrib`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)
