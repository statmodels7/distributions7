# Skew Normal Second Response Derivative in the Centered Parametrization

Computes \\\partial^2\ell/\partial y^2\\ by delegating to
[`distrib_hess_y.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.SkewNormal1Distrib.md)
at the implied direct parameters. Like the first response derivative it
is the parent's unchanged, and is defined at \\\gamma_1 = 0\\.

The value is strictly negative at every observation and every skewness:
the skew normal log-density is concave in the response.

## Arguments

- distrib:

  A `SkewNormal2Distrib` object, from
  [`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma` and `gamma1`, each a
  numeric vector of length 1 or of the length of `y`.

- ...:

  Passed to
  [`distrib_hess_y.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.SkewNormal1Distrib.md).

## Value

A numeric vector of the length of the recycled inputs, negative
throughout.

## See also

[`distrib_grad_y.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.SkewNormal2Distrib.md)
for the first derivative,
[`distrib_hess_y.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.SkewNormal1Distrib.md)
for the closed form, and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- skewnormal2_distrib()
y <- c(-1, 0.3, 1.7)
th <- list(mu = 0, sigma = 1, gamma1 = 0.5)

# Against a central difference of the response derivative.
eps <- 1e-5
rbind(analytic = distrib_hess_y(d, y, th),
      numeric = (distrib_grad_y(d, y + eps, th) -
                 distrib_grad_y(d, y - eps, th)) / (2 * eps))
#>               [,1]       [,2]       [,3]
#> analytic -1.863173 -0.7190311 -0.4753314
#> numeric  -1.863173 -0.7190311 -0.4753314

# Concave in the response at every skewness the family reaches.
vapply(c(-0.9, -0.3, 0.3, 0.9), function(g)
  max(distrib_hess_y(d, seq(-6, 6, by = 0.5),
                     list(mu = 0, sigma = 1, gamma1 = g))), 0)
#> [1] -0.3790360 -0.5594095 -0.5594095 -0.3790360
```
