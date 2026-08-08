# Laplace Mixed Derivatives in Location and Rate

Closed form, away from the kink at \\y = \mu\\: with \\\ell^{(y)} =
-\lambda\\\mathrm{sign}(y-\mu)\\, the location component is 0 almost
everywhere and the rate component is \\-\mathrm{sign}(y-\mu)\\.

## Arguments

- distrib:

  A `Laplace2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `lambda`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with one numeric vector per parameter.
