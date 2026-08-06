# CDF Derivatives of a Weibull

Closed form at both orders. With \\u = (q/\mu)^\sigma\\, \\L =
\log(q/\mu)\\ and \\F = 1 - e^{-u}\\, every derivative of \\F\\ follows
from those of \\u\\, which are elementary: \\\partial u/\partial\mu =
-\sigma u/\mu\\ and \\\partial u/\partial\sigma = uL\\.

## Usage

``` r
weibull_cdf_deriv(distrib, q, theta, order)
```

## Arguments

- distrib:

  A `WeibullDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu` and `sigma`.

- order:

  The derivative order, 1 or 2.

## Value

A named list of derivative components of \\F\\.

## See also

[`weibull1_distrib`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)
