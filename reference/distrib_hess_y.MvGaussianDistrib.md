# Multivariate Gaussian Response Hessian

Computes the second derivative of the log-density in the response,
\$\$\frac{\partial^2 \ell}{\partial y\\ \partial y^\top} =
-\Sigma^{-1}.\$\$ The log-density is quadratic in \\y\\, so the matrix
is the same at every observation and the response is not read. One \\p
\times p\\ matrix is returned for the whole sample, and not one matrix
per row.

## Arguments

- distrib:

  An
  [MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object, from
  [`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md).

- y:

  An \\n \times p\\ numeric matrix of observations. It is not read.

- theta:

  A named list of parameters in which each component is a single number.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A \\p \times p\\ numeric matrix, symmetric and negative definite.

## See also

[`distrib_grad_y.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.MvGaussianDistrib.md)
for the first derivative,
[`distrib_cross2_y.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.MvGaussianDistrib.md)
for its derivative in the parameters, and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
set.seed(1)
y <- distrib_rng(d, 4, theta)

distrib_hess_y(d, y, theta)
#>            [,1]       [,2]
#> [1,] -1.0141552  0.5399435
#> [2,]  0.5399435 -1.4918247
-solve(mv_sigma(d, theta))
#>            v1         v2
#> v1 -1.0141552  0.5399435
#> v2  0.5399435 -1.4918247

# One matrix is returned whatever the sample size, and it is the same at
# every observation.
identical(distrib_hess_y(d, y, theta), distrib_hess_y(d, y[1, ], theta))
#> [1] TRUE

# The result is compared with a numerical Hessian at one observation.
max(abs(distrib_hess_y(d, y, theta) -
        numDeriv::hessian(function(z) distrib_pdf(d, z, theta, log = TRUE),
                          y[1, ])))
#> [1] 2.747602e-11
```
