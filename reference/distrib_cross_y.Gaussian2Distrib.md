# Gaussian Mixed Derivatives in Mean and Variance

Closed form: \\\ell^{(y)} = -(y-\mu)/\sigma^2\\ gives \\1/\sigma^2\\ and
\\(y-\mu)/\sigma^4\\.

## Arguments

- distrib:

  A `Gaussian2Distrib` object.

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
