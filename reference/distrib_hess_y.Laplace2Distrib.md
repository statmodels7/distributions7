# Laplace Second Derivative in the Response, Rate Parametrization

Returns zero for every observation. The log-density \\\log(\lambda/2) -
\lambda\|y-\mu\|\\ is linear in \\y\\ on each side of the location, so
its second derivative in the response vanishes wherever it exists. At
\\y = \mu\\ the first derivative drops from \\\lambda\\ to \\-\lambda\\
and the second derivative does not exist; the returned zero is the value
away from that single point.

## Arguments

- distrib:

  A `Laplace2Distrib` object, from
  [`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with components `mu` and `lambda`. Neither is read.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of zeros, of length `length(y)`.

## See also

[`distrib_grad_y.Laplace2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Laplace2Distrib.md)
for the first derivative, which drops at the location;
[`distrib_hess_y.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.LaplaceDistrib.md)
for the scale parametrization;
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- laplace2_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, lambda = 2)

distrib_hess_y(d, y, th)
#> [1] 0 0 0

# A location family: the same as the curvature in the location.
all.equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
#> [1] TRUE

# The first derivative drops by 2 lambda across the location.
diff(distrib_grad_y(d, 0.4 + c(-1e-9, 1e-9), th))
#> [1] -4
-2 * 2
#> [1] -4
```
