# Skew Normal Cumulative Distribution Function

Computes the distribution function through Owen's T function: \$\$F(q;
\mu, \sigma, \alpha) = \Phi(z) - 2\\T(z, \alpha), \qquad z = (q -
\mu)/\sigma.\$\$

## Arguments

- distrib:

  A `SkewNormal1Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu`, `sigma` and `alpha`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le q)\\,
  otherwise \\P(Y \> q)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## Value

A numeric vector of cumulative probabilities.

## Details

The identity is Azzalini's. Evaluating it costs one bounded
one-dimensional quadrature per observation, which is cheaper and more
accurate than the base class's route of integrating the density over a
semi-infinite range.

## See also

[`skewnormal1_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
