# Multivariate Gaussian Fourth Derivative in Two Responses and Two Parameters

Computes \\\partial^4\ell/\partial y\\\partial y^\top
\partial\theta_a\partial\theta_b\\, one \\p \times p\\ matrix per
unordered pair of parameters. The response Hessian \\-\Sigma^{-1}\\ does
not involve the mean, so a pair naming any mean parameter gives the zero
matrix; for a pair of free values, \$\$\frac{\partial^4\ell}{\partial
y\\\partial y^\top \partial\eta_k\partial\eta_l} =
\Sigma^{-1}A\_{kl}\Sigma^{-1} - \Sigma^{-1}\\\left(A_l\Sigma^{-1}A_k +
A_k\Sigma^{-1}A_l\right)\\\Sigma^{-1},\$\$ with \\A_k\\ and \\A\_{kl}\\
the first and second derivatives of \\\Sigma\\. No component depends on
the observation, so `y` is not read.

## Arguments

- distrib:

  An
  [MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object, from
  [`mvgaussian_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian_distrib.md).

- y:

  An \\n \times p\\ numeric matrix of observations. Not read.

- theta:

  A named list of parameters, each component a single number.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch. The two coincide here.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of \\p \times p\\ numeric matrices, keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md).
Every key naming a mean parameter holds the zero matrix.

## Notation

\\\Sigma\\ is the covariance, \\\eta\\ the free vector of the matrix
parametrization, \\A_k\\ and \\A\_{kl}\\ its first and second derivative
arrays, and \\\ell\\ the log-density of one observation.

## See also

[`distrib_cross2_y.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.MvGaussianDistrib.md)
for the same derivative one parameter down,
[`distrib_grad_y_hess.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.MvGaussianDistrib.md)
for its sibling with one response index, and
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md)
for the generic.

## Examples

``` r
d <- mvgaussian_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
y <- matrix(0, 4, 2)

hh <- distrib_hess_y_hess(d, y, theta)
length(hh) == length(hess_names(d@params))
#> [1] TRUE

# Any pair naming a mean gives the zero matrix.
hh$mu1_sigma_L2.1
#>      [,1] [,2]
#> [1,]    0    0
#> [2,]    0    0

# A matrix pair against a second difference of the response Hessian.
h <- 1e-4
f <- function(a, b) {
  t2 <- theta
  t2$sigma_log_L1 <- t2$sigma_log_L1 + a
  t2$sigma_L2.1 <- t2$sigma_L2.1 + b
  distrib_hess_y(d, y, t2)
}
num <- (f(h, h) - f(h, -h) - f(-h, h) + f(-h, -h)) / (4 * h * h)
round(hh$sigma_log_L1_sigma_L2.1, 6)
#>           [,1]      [,2]
#> [1,]  1.954244 -1.349859
#> [2,] -1.349859  0.000000
round(num, 6)
#>           [,1]      [,2]
#> [1,]  1.954244 -1.349859
#> [2,] -1.349859  0.000000
```
