# NB1 Analytical Expected Hessian

Every term carrying \\P\\ drops out, its expectation vanishing by the
first Bartlett identity, and what remains needs only
\\\mathbb{E}\[\psi'(Y+r)\]\\. That has no closed form and is summed
against the exact mass to a far-tail quantile, as the NB2 kernel does;
`approx` is therefore ignored.

## Arguments

- distrib:

  A `NegBin1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `theta`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  Ignored.

- nsim:

  Ignored.

- ...:

  Unused.

- threads:

  How many threads the kernel may use; below the measured internal
  threshold it stays sequential whatever the count says.

## Value

A named list of expected second-derivative components.

## See also

[`negbin1_distrib`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md)
