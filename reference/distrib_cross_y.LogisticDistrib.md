# Logistic Mixed Derivatives

Closed form from the location-scale identity: the location component is
\\-\ell^{(yy)}\\ and the scale one \\-z\ell^{(yy)} -
\ell^{(y)}/\sigma\\.

## Arguments

- distrib:

  A `LogisticDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `sigma`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with one numeric vector per parameter.
