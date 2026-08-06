# Gamma Mixed Derivatives in Mean and Dispersion

Closed form. The shape is \\1/\phi\\ and the rate \\1/(\phi\mu)\\, so
the mean component is \\1/(\phi\mu^2)\\ and the dispersion one
\\1/(\phi^2\mu) - 1/(\phi^2 y)\\.

## Arguments

- distrib:

  A `Gamma1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `phi`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with one numeric vector per parameter.
