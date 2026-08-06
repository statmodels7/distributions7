# Lognormal Mixed Derivatives

Closed form. On the log scale the family is location-scale, and
\\\ell^{(y)} = -1/y - (\log y - \mu)/(\sigma^2 y)\\ gives \\1/(\sigma^2
y)\\ and \\(\log y - \mu)/(\sigma^4 y)\\.

## Arguments

- distrib:

  A `Lognormal1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `sigma2`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with one numeric vector per parameter.
