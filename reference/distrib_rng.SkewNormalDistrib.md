# Skew Normal Random Number Generator

Generates draws exactly, from the stochastic representation \\Z =
\delta\|U_0\| + \sqrt{1-\delta^2}\\U_1\\ with \\U_0, U_1\\ independent
standard normal and \\\delta = \alpha/\sqrt{1+\alpha^2}\\. No inversion
or rejection is involved.

## Arguments

- distrib:

  A `SkewNormalDistrib` object.

- n:

  Number of observations to generate.

- theta:

  A list containing `mu`, `sigma` and `alpha`.

## Value

A numeric vector of random draws.

## See also

[`skewnormal_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal_distrib.md)
