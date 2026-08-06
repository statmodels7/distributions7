# Generalized Gamma Mixed Derivatives

Closed form. With \\w = (y/a)^p\\ and \\L = \log(y/a)\\, \\\ell^{(y)} =
((d-1) - pw)/y\\, so the components are \\p^2w/(ay)\\, \\1/y\\ and
\\-w(1 + pL)/y\\.

## Arguments

- distrib:

  A `GenGamma1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `a`, `d` and `p`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with one numeric vector per parameter.
