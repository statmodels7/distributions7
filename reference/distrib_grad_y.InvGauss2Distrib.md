# Inverse Gaussian Response and Mixed Derivatives in Mean and Shape

Obtains the derivatives from those of the dispersion parametrization, at
\\\phi = 1/\lambda\\. A derivative in the response is unaffected by the
parametrization, so the response derivatives are the parent's unchanged,
and the mixed block is one first-order chain.

## Arguments

- distrib:

  An `InvGauss2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `lambda`.

- ...:

  Unused.

- scale:

  Handled by the generic after dispatch; this method returns the
  parameter scale.

## Value

A numeric vector.
