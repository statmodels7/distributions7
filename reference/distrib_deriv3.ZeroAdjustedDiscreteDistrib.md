# Hurdle Third Derivatives

The likelihood separates, so mixed components vanish at every order; the
\\\theta\\ part is the parent's derivative less that of the truncation
constant \\\log(1-f_0)\\.

## Arguments

- distrib:

  A `ZeroAdjustedDiscreteDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with the parent's parameters followed by `za`.

- expected:

  Logical; if `TRUE`, the expected derivatives.

## Value

A named list of derivative components.

## See also

[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
