# Skew t Fifth Derivatives

Computes the fifty-six fifth derivatives of the log-density with the
compiled kernel described on
[`distrib_gradient.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewTDistrib.md).

## Arguments

- distrib:

  A `SkewTDistrib` object, from
  [`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma`, `alpha` and `nu`.

- scale:

  Either `"parameter"`, the default, or `"link"`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use.

## Value

A named list of fifty-six numeric vectors, one per distinct fifth-order
component, as
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
names them.

## See also

[`distrib_deriv4.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.SkewTDistrib.md)
for the order below and
[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md)
for the generic.

## Examples

``` r
d <- skewt_distrib()
d5 <- distrib_deriv5(d, c(-1.5, 0.4, 2.1),
                     list(mu = 0, sigma = 1, alpha = 3, nu = 6))
length(d5)
#> [1] 56
```
