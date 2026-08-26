# Skew Normal Response Derivative

Computes \\\partial\ell/\partial y\\, the derivative of the log-density
in the response. With \\z = (y-\mu)/\sigma\\ and \\R\\ the inverse Mills
ratio at \\t = \alpha z\\, \$\$\dfrac{\partial \ell}{\partial y} =
\dfrac{\alpha R - z}{\sigma}.\$\$ This is minus the derivative in
\\\mu\\, exactly, because the response and the location enter the
density only through their difference. The identity holds for every
location family, so this method is one sign away from
[`distrib_gradient.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewNormal1Distrib.md).

The quantity is what a censored likelihood and a quantile residual need,
and it is the first ingredient of
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md).

## Arguments

- distrib:

  A `SkewNormal1Distrib` object, from
  [`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma` and `alpha`, each a numeric
  vector of length 1 or of the length of `y`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of the length of the recycled inputs, one value per
observation.

## Notation

\\\ell\\ is the log-density of one observation, \\z = (y-\mu)/\sigma\\
and \\R(t) = \phi(t)/\Phi(t)\\ the inverse Mills ratio.

## See also

[`distrib_hess_y.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.SkewNormal1Distrib.md)
for the second derivative,
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the mixed one,
[`distrib_gradient.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewNormal1Distrib.md)
for the parameter derivatives, and
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic.

## Examples

``` r
d <- skewnormal1_distrib()
y <- c(-1.5, -0.3, 0.4, 2.1)
th <- list(mu = 0, sigma = 1, alpha = 3)

# It is minus the score in the location, exactly.
all.equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
#> [1] TRUE

# Against a central difference of the log-density.
eps <- 1e-6
rbind(analytic = distrib_grad_y(d, y, th),
      numeric = (distrib_pdf(d, y + eps, th, log = TRUE) -
                 distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps))
#>              [,1]     [,2]      [,3] [,4]
#> analytic 15.61296 4.636929 0.2583096 -2.1
#> numeric  15.61296 4.636929 0.2583096 -2.1
```
