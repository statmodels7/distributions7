# Skew Normal Second Response Derivative

Computes \\\partial^2\ell/\partial y^2\\, the curvature of the
log-density in the response. With \\R' = -R(t+R)\\ at \\t = \alpha z\\,
\$\$\dfrac{\partial^2 \ell}{\partial y^2} = \dfrac{\alpha^2 R' -
1}{\sigma^2},\$\$ which is the same expression as
\\\partial^2\ell/\partial\mu^2\\: two signs cancel where one did not at
first order, so the response curvature equals the location curvature
instead of being its negative.

Since \\R'\\ is negative for every \\t\\, the value is strictly negative
at every observation and every shape. The skew normal log-density is
concave in the response, as the Gaussian's is.

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

A numeric vector of the length of the recycled inputs, negative
throughout.

## Notation

\\z = (y-\mu)/\sigma\\, \\t = \alpha z\\, \\R\\ the inverse Mills ratio
and \\R' = -R(t+R)\\ its derivative.

## See also

[`distrib_grad_y.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.SkewNormal1Distrib.md)
for the first derivative,
[`distrib_hessian.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.SkewNormal1Distrib.md)
for the parameter curvature it shares an expression with, and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- skewnormal1_distrib()
y <- c(-1.5, -0.3, 0.4, 2.1)
th <- list(mu = 0, sigma = 1, alpha = 3)

# It equals the curvature in the location, without a sign change.
all.equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
#> [1] TRUE

# Against a central difference of the response derivative.
eps <- 1e-5
rbind(analytic = distrib_hess_y(d, y, th),
      numeric = (distrib_grad_y(d, y + eps, th) -
                 distrib_grad_y(d, y - eps, th)) / (2 * eps))
#>               [,1]      [,2]      [,3] [,4]
#> analytic -9.650673 -8.099247 -3.803286   -1
#> numeric  -9.650673 -8.099247 -3.803286   -1

# Concave in the response at every shape, and bounded above by
# -1/sigma^2, which is the Gaussian's own curvature.
vapply(c(-20, -1, 0, 1, 20), function(a)
  max(distrib_hess_y(d, seq(-6, 6, by = 0.5), list(mu = 0, sigma = 1, alpha = a))), 0)
#> [1] -1 -1 -1 -1 -1
```
