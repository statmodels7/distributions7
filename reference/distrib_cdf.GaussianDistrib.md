# Gaussian Cumulative Distribution Function

Computes the cumulative distribution function for the Gaussian
distribution: \$\$F(q; \mu, \sigma) =
\Phi\left(\dfrac{q-\mu}{\sigma}\right)\$\$

## Arguments

- distrib:

  A `GaussianDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing the parameters `mu` and `sigma`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le q)\\,
  otherwise \\P(Y \> q)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## See also

[`gaussian_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian_distrib.md)
