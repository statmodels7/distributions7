# Response Hessian of a Multivariate Family Without a Method

Signals an error, for the reason that
[`distrib_grad_y.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.multivariate_distrib.md)
gives: the univariate numerical fallback differences along a line and
would produce a scalar where a \\p \times p\\ matrix is wanted. A family
that can supply \\\partial^2\ell/\partial y\\\partial y^\top\\ registers
its own method.

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object with no method of its own.

- y:

  An \\n \times p\\ numeric matrix of observations. It is not read.

- theta:

  A named list of parameters. It is not read.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

Does not return: it signals an error naming the family.

## Details

The two families with a method return different shapes: the gaussian's
Hessian does not depend on the observation and is returned as one \\p
\times p\\ matrix, while the Student t's does and is returned as a \\p
\times p \times n\\ array.

## See also

[`distrib_hess_y.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.MvGaussianDistrib.md)
and
[`distrib_hess_y.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.MvStudentTDistrib.md)
for the two families with a method, and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- dirichlet_distrib(3)
theta <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8)
set.seed(1)
y <- distrib_rng(d, 3, theta)
try(distrib_hess_y(d, y, theta))
#> Error : distrib_hess_y() is not defined for 'dirichlet [3d, mean=simplex]': the univariate numerical fallback differences along a line; register a closed form on the family.

# The two families with a method return different shapes.
th <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
           sigma_L2.1 = 0.5)
yy <- rbind(c(1, -1), c(0, 0))
dim(distrib_hess_y(mvgaussian1_distrib(2), yy, th))
#> [1] 2 2
dim(distrib_hess_y(mvstudent_t1_distrib(2), yy, c(th, list(nu = 6))))
#> [1] 2 2 2
```
