# Skew Normal Log-CDF Derivatives

Closed form at every order, in the shape as well as in the location and
the scale. With \\z = (q-\mu)/\sigma\\ the distribution function is
\\\Phi(z) - 2T(z, \alpha)\\, and Owen's \\T\\ has elementary partial
derivatives, so the integral in its definition is differentiated away at
the first order and never has to be differentiated again. The location
and the scale then enter only through \\z\\, by the same chain rule the
other location-scale families use.

## Usage

``` r
sn_cdf_deriv_k(distrib, q, theta, order)
```

## Arguments

- distrib:

  A `SkewNormal1Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu`, `sigma` and `alpha`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

- ...:

  Unused.

## Value

A named list, one vector per component.

## See also

[`skewnormal1_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
