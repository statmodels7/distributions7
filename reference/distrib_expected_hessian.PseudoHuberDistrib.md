# Pseudo-Huber Analytical Expected Hessian

Computes the expected Hessian of the Pseudo-Huber log-density. No closed
form exists, so the non-zero components are evaluated by numerical
integration of the observed Hessian against the density (via
[`expectation`](https://statmodels7.github.io/distributions7/reference/expectation.md)).
By symmetry \\\mathbb{E}\[r\] = \mathbb{E}\[r^3\] = 0\\, hence the
\\\mu\sigma\\ and \\\mu\nu\\ components are exactly 0.

## Arguments

- distrib:

  A `PseudoHuberDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu`, `sigma` and `nu`.

## Value

A list containing the vectors of expected second derivatives.

## See also

[`pseudohuber_distrib`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
