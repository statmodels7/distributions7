# Beta-Binomial Mass Function in Its Shapes

\$\$P(Y = y) = \binom{n}{y}\dfrac{B(y+\alpha,
n-y+\beta)}{B(\alpha,\beta)}\$\$

## Arguments

- distrib:

  A `BetaBinom2Distrib` object.

- y:

  A numeric vector of counts.

- theta:

  A list with `alpha` and `beta`.

- log:

  Logical; if `TRUE`, returns the log-probability.

## Value

A numeric vector.

## See also

[`betabinom2_distrib`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md)
