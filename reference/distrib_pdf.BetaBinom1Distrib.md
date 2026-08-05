# Beta-Binomial Probability Mass Function

\$\$P(Y = y) = \binom{n}{y} \dfrac{B(y + \alpha,\\ n - y +
\beta)}{B(\alpha, \beta)}, \qquad \alpha = \dfrac{\mu}{\sigma}, \quad
\beta = \dfrac{1-\mu}{\sigma}\$\$

## Arguments

- distrib:

  A `BetaBinom1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `sigma`.

- log:

  Logical; if `TRUE`, returns the log-probability.

## Value

A numeric vector of probability values.

## See also

[`betabinom1_distrib`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md)
