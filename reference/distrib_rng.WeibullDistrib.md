# Weibull Random Number Generator

Generates random numbers by inverse transform, which is exact here
because the quantile function is elementary.

## Arguments

- distrib:

  A `WeibullDistrib` object.

- n:

  Number of observations to generate.

- theta:

  A list containing the parameters `mu` and `sigma`.

## Value

A numeric vector of random draws.

## See also

[`weibull_distrib`](https://statmodels7.github.io/distributions7/reference/weibull_distrib.md)
