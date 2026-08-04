# Zero-Inflated Probability Mass Function

\$\$P(Y=y) = \zeta\\\mathbb{I}(y=0) + (1-\zeta) f(y; \theta)\$\$ where
\\f\\ is the PMF of the wrapped distribution.

## Arguments

- distrib:

  A `ZeroInflatedDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with the parent's parameters followed by `zi`.

- log:

  Logical; if `TRUE`, returns the log-probability.

## Value

A numeric vector of density values.

## See also

[`zero_inflated`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
