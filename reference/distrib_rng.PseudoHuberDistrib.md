# Pseudo-Huber Random Number Generator

Generates random numbers from the Pseudo-Huber distribution via inverse
transform sampling (numerical quantile function applied to uniform
draws).

## Arguments

- distrib:

  A `PseudoHuberDistrib` object.

- n:

  Number of observations to generate.

- theta:

  A list containing the parameters `mu`, `sigma` and `nu`.

## Value

A numeric vector of random draws.

## See also

[`pseudohuber_distrib`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
