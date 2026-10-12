# Response Gradient of a Multivariate Family Without a Method

Signals an error. The univariate numerical fallback differences the
log-density along a line, which in \\p\\ dimensions gives a directional
derivative, a number of the wrong shape. A family that can supply
\\\partial\ell/\partial y\\ registers its own method, as the gaussian
and the Student t do.

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object with no method of its own.

- y:

  An \\n \times p\\ numeric matrix of observations. Not read.

- theta:

  A named list of parameters. Not read.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

Does not return: it signals an error naming the family. A family that
registers a method returns an \\n \times p\\ matrix.

## Details

[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
tests for a method of the family's own with
[`has_mv_grad_y()`](https://statmodels7.github.io/distributions7/reference/has_mv_grad_y.md)
and leaves out the response-derivative check when there is none.

## See also

[`distrib_grad_y.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.MvGaussianDistrib.md)
and
[`distrib_grad_y.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.MvStudentTDistrib.md)
for the two families with a method,
[`distrib_hess_y.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.multivariate_distrib.md)
for the second order, and
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic.

## Examples

``` r
# The Dirichlet registers no response gradient, so the method of the class
# signals an error.
d <- dirichlet_distrib(3)
theta <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8)
set.seed(1)
y <- distrib_rng(d, 3, theta)
try(distrib_grad_y(d, y, theta))
#> Error : distrib_grad_y() is not defined for 'dirichlet [3d, mean=simplex]': the univariate numerical fallback differences along a line; register a closed form on the family.

# The gaussian registers one, which returns an n by p matrix.
g <- mvgaussian1_distrib(2)
th <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
           sigma_L2.1 = 0.5)
dim(distrib_grad_y(g, rbind(c(1, -1), c(0, 0)), th))
#> [1] 2 2
```
