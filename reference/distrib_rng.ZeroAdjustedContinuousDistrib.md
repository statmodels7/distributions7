# Zero-Adjusted Continuous Random Number Generator

Draws zeros with probability \\\pi\\ and otherwise samples from the
parent distribution.

## Arguments

- distrib:

  A `ZeroAdjustedContinuousDistrib` object.

- n:

  Number of observations to generate.

- theta:

  A list with the parent's parameters followed by `za`.

## Value

A numeric vector of random draws.

## See also

[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
