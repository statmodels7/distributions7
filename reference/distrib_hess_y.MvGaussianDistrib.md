# Multivariate Gaussian Response Hessian

Computes the second derivative of the log-density in the response,
\$\$\frac{\partial^2 \ell}{\partial y\\ \partial y^\top} =
-\Sigma^{-1}.\$\$ The quadratic form is quadratic in \\y\\, so the
matrix is the same at every observation and the response is read only to
confirm its shape. One \\p \times p\\ matrix is returned, not \\n\\
copies of it, and a consumer that expects one matrix per row must handle
this family's shape explicitly.

## Arguments

- distrib:

  An
  [MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object, from
  [`mvgaussian_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian_distrib.md).

- y:

  An \\n \times p\\ numeric matrix of observations. Its values do not
  enter the result.

- theta:

  A named list of parameters, each component a single number.

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
d <- mvgaussian_distrib(2)
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

# One matrix whatever the sample size, and it does not move with y.
identical(distrib_hess_y(d, y, theta), distrib_hess_y(d, y[1, ], theta))
#> [1] TRUE

# Against a numerical Hessian at one observation.
max(abs(distrib_hess_y(d, y, theta) -
        numDeriv::hessian(function(z) distrib_pdf(d, z, theta, log = TRUE),
                          y[1, ])))
#> [1] 2.747602e-11
```
