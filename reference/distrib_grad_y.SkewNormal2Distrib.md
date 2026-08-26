# Skew Normal Response Derivative in the Centered Parametrization

Computes \\\partial\ell/\partial y\\ by delegating to
[`distrib_grad_y.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.SkewNormal1Distrib.md)
at the implied direct parameters. The value is the parent's unchanged: a
change of parameters does not touch a derivative in the response, the
two variables being separate arguments of the same log-density.

It is therefore defined at \\\gamma_1 = 0\\, where the parameter
derivatives are not: nothing here differentiates the map.

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
  [`distrib_grad_y.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.SkewNormal1Distrib.md).

## Value

A numeric vector of the length of the recycled inputs.

## See also

[`distrib_hess_y.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.SkewNormal2Distrib.md)
for the second derivative,
[`distrib_grad_y.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.SkewNormal1Distrib.md)
for the closed form it delegates to, and
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic.

## Examples

``` r
d <- skewnormal2_distrib()
y <- c(-1, 0.3, 1.7)
th <- list(mu = 0, sigma = 1, gamma1 = 0.5)

# Against a central difference of the log-density in the response.
eps <- 1e-6
rbind(analytic = distrib_grad_y(d, y, th),
      numeric = (distrib_pdf(d, y + eps, th, log = TRUE) -
                 distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps))
#>              [,1]       [,2]     [,3]
#> analytic 1.096524 -0.5631452 -1.30601
#> numeric  1.096524 -0.5631452 -1.30601

# Defined at zero skewness, where the parameter derivatives are not.
all.equal(distrib_grad_y(d, y, list(mu = 0, sigma = 1, gamma1 = 0)), -y)
#> [1] TRUE
```
