# Elastic-Net Second Response Derivative

Computes \\\partial^2\ell/\partial y^2 = -c\\, with \\c =
\lambda(1-\alpha)\\: the absolute value contributes nothing away from
the location, so the curvature in the response is the Gaussian part
alone. It is constant, free of the data, and strictly negative.

At the location itself the second derivative carries a point mass
\\-2a\delta(y-\mu)\\, which this value omits, exactly as
[`distrib_hessian.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.EnetDistrib.md)'s
`mu_mu` does.

Every response derivative of **third order and beyond is zero**: the
log-density is quadratic in \\y\\ away from the kink.

## Arguments

- distrib:

  An `EnetDistrib` object, from
  [`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md).

- y:

  A numeric vector of observations. Its values do not enter the result;
  only its length does.

- theta:

  A named list with components `mu`, `lambda` and `alpha`, each a
  numeric vector of length 1 or of the length of `y`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of the length of `y`, every entry \\-c\\.

## Notation

\\c = \lambda(1-\alpha)\\ is the Gaussian rate and \\a = \lambda\alpha\\
the Laplace one.

## See also

[`distrib_grad_y.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.EnetDistrib.md)
for the first derivative,
[`distrib_hessian.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.EnetDistrib.md)
for the parameter curvature it shares an expression with, and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- enet_distrib()
y <- c(-1.5, -0.3, 0.4, 2.1)
th <- list(mu = 0, lambda = 2, alpha = 0.5)

# Constant, and equal to the curvature in the location.
c(hess_y = unique(distrib_hess_y(d, y, th)),
  mu_mu = unique(distrib_hessian(d, y, th)$mu_mu),
  minus_c = -(2 * (1 - 0.5)))
#>  hess_y   mu_mu minus_c 
#>      -1      -1      -1 

# Against a central difference of the response derivative, away from the
# kink.
eps <- 1e-5
rbind(analytic = distrib_hess_y(d, y, th),
      numeric = (distrib_grad_y(d, y + eps, th) -
                 distrib_grad_y(d, y - eps, th)) / (2 * eps))
#>          [,1] [,2] [,3] [,4]
#> analytic   -1   -1   -1   -1
#> numeric    -1   -1   -1   -1

# The Gaussian rate alone, so it goes to zero as alpha approaches one.
vapply(c(0.1, 0.5, 0.9, 0.999),
       function(a) distrib_hess_y(d, 1, list(mu = 0, lambda = 2,
                                             alpha = a)), 0)
#> [1] -1.800 -1.000 -0.200 -0.002
```
