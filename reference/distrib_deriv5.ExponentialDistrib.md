# Exponential Fifth-Order Derivatives

Computes the fifth derivative of the exponential log-density \\\ell =
-\log\mu - y/\mu\\ with respect to its mean, in closed form from a
compiled kernel, \$\$\ell^{(\mu\mu\mu\mu\mu)} = \dfrac{120 y}{\mu^6} -
\dfrac{24}{\mu^5}.\$\$ With `scale = "link"` the components are composed
with the link derivatives by the routine that every family with a
fifth-order kernel uses.

## Arguments

- distrib:

  An `ExponentialDistrib` object, from
  [`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of the length of `y`.

- scale:

  One of `"parameter"` (the default) or `"link"`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use.
  Defaults to `1L`.

## Value

A named list with one numeric vector, `mu_mu_mu_mu_mu`, of length
`max(length(y), length(mu))`.

## See also

[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the order below and
[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md)
for the generic.

## Examples

``` r
d <- exponential_distrib()
y <- c(0.2, 1, 4)
all.equal(distrib_deriv5(d, y, list(mu = 1.5))$mu_mu_mu_mu_mu,
          120 * y / 1.5^6 - 24 / 1.5^5)
#> [1] TRUE
```
