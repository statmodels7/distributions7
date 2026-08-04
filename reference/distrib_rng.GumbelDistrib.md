# Gumbel Random Number Generator

Generates random numbers as \\\mu - \sigma \log E\\ with \\E\\ standard
exponential, which is inverse transform written without a logarithm of a
uniform.

## Arguments

- distrib:

  A `GumbelDistrib` object.

- n:

  Number of observations to generate.

- theta:

  A list containing the parameters `mu` and `sigma`.

## Value

A numeric vector of random draws.

## See also

[`gumbel_distrib`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md)
