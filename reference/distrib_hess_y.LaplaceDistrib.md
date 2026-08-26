# Laplace Second Derivative in the Response

Returns zero for every observation. The Laplace log-density is
\\-\log(2\sigma) - \|y-\mu\|/\sigma\\, which is **linear** in \\y\\ on
each side of the location, so its second derivative in the response
vanishes wherever it exists. At \\y = \mu\\ the first derivative drops
from \\1/\sigma\\ to \\-1/\sigma\\ and the second derivative does not
exist; the returned zero is the value away from that single point.

The same fact appears in the parameters as the zero `mu_mu` component of
[`distrib_hessian.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.LaplaceDistrib.md),
the family being a location family.

## Arguments

- distrib:

  A `LaplaceDistrib` object, from
  [`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with components `mu` and `sigma`. Neither is read.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of zeros, of length `length(y)`.

## See also

[`distrib_grad_y.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.LaplaceDistrib.md)
for the first derivative, which jumps at the location;
[`distrib_hessian.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.LaplaceDistrib.md)
for the same vanishing curvature in the parameters;
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- laplace_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)

distrib_hess_y(d, y, th)
#> [1] 0 0 0

# A location family: the same as the curvature in the location.
all.equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
#> [1] TRUE

# The first derivative drops by 2/sigma across the location, which is the
# curvature the zero above does not see.
diff(distrib_grad_y(d, 0.4 + c(-1e-9, 1e-9), th))
#> [1] -1.333333
-2 / 1.5
#> [1] -1.333333
```
