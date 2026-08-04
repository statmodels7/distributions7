# Inverse-Gaussian Cumulative Distribution Function

Computes the cumulative distribution function for the Inverse-Gaussian
distribution: \$\$F(q; \mu, \phi) = \Phi\\\left(\sqrt{\dfrac{1}{\phi
q}}\left(\dfrac{q}{\mu}-1\right)\right) + e^{2/(\phi\mu)}\\
\Phi\\\left(-\sqrt{\dfrac{1}{\phi
q}}\left(\dfrac{q}{\mu}+1\right)\right)\$\$ where \\\Phi\\ is the
standard normal CDF.

## Arguments

- distrib:

  An `InvGaussDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing the parameters `mu` and `phi`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le q)\\,
  otherwise \\P(Y \> q)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## Value

A numeric vector of cumulative probabilities.

## See also

[`invgauss_distrib`](https://statmodels7.github.io/distributions7/reference/invgauss_distrib.md)
