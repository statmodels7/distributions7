# Skew t Random Number Generator

Generates draws exactly, from the scale mixture \\Y = \mu + \sigma
Z/\sqrt{V/\nu}\\ with \\Z\\ standard skew normal of shape \\\alpha\\ and
\\V \sim \chi^2\_\nu\\ independent of it. No inversion or rejection is
involved.

## Arguments

- distrib:

  A `SkewTDistrib` object.

- n:

  Number of observations to generate.

- theta:

  A list containing `mu`, `sigma`, `alpha` and `nu`.

## Value

A numeric vector of random draws.

## See also

[`skewt_distrib`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md)
