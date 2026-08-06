# Gaussian Mixed Derivatives in Mean and Precision

Closed form: \\\ell^{(y)} = -\tau(y-\mu)\\ gives \\\tau\\ and
\\-(y-\mu)\\.

## Arguments

- distrib:

  A `Gaussian3Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `tau`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with one numeric vector per parameter.
