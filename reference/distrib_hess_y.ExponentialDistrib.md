# Exponential Second Derivative in the Response

Returns zero for every observation. The exponential log-density is
linear in \\y\\ on the whole support, so its second derivative in the
response vanishes everywhere; there is no kink to qualify the statement,
unlike the Laplace, whose zero holds only away from a point.

## Arguments

- distrib:

  An `ExponentialDistrib` object, from
  [`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with the single component `mu`, which is not read.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of zeros, of length `length(y)`.

## See also

[`distrib_grad_y.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.ExponentialDistrib.md)
for the constant first derivative;
[`distrib_hess_y.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.LaplaceDistrib.md),
which is zero for a different reason;
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- exponential_distrib()
th <- list(mu = 2)

distrib_hess_y(d, c(0.3, 1.1, 4.0), th)
#> [1] 0 0 0

# The first derivative is constant, so its difference is zero everywhere,
# with no exceptional point.
eps <- 1e-6
y <- c(1e-8, 0.5, 5, 50)
(distrib_grad_y(d, y + eps, th) - distrib_grad_y(d, y - eps, th)) / (2 * eps)
#> [1] 0 0 0 0
```
