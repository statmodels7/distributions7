# Multivariate Gaussian Third and Fourth Response Derivatives

Returns the third and fourth derivatives of the log-density in the
response, and the derivative of the third one in each parameter. The
log-density is quadratic in \\y\\, with response Hessian
\\-\Sigma^{-1}\\, so all three quantities are zero. The arrays follow
[`distrib_hess_y.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.MvGaussianDistrib.md):
they are constant in \\y\\ and carry no observation dimension.

## Arguments

- distrib:

  An
  [MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object, from
  [`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md)
  or
  [`mvgaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md).

- y:

  An \\n \times p\\ numeric matrix of observations. A vector of length
  \\p\\ is read as a single observation.

- theta:

  A named list of parameters in which each component is a single number.

- scale:

  One of `"parameter"` (the default) or `"link"`. Every link of this
  family is the identity, so the two scales coincide.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md)
returns a \\p \times p \times p\\ array of zeros,
[`distrib_deriv4_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md)
a \\p \times p \times p \times p\\ one, and
[`distrib_cross3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross3_y.md)
a named list of \\p \times p \times p\\ arrays of zeros, one per
parameter, in `distrib@params` order.

## See also

[`distrib_hess_y.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.MvGaussianDistrib.md)
for the order below, and
[`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md)
and
[`distrib_cross3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross3_y.md)
for the generics.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
set.seed(1)
y <- distrib_rng(d, 4, theta)
dim(distrib_deriv3_y(d, y, theta))
#> [1] 2 2 2
dim(distrib_deriv4_y(d, y, theta))
#> [1] 2 2 2 2
names(distrib_cross3_y(d, y, theta))
#> [1] "mu1"          "mu2"          "sigma_log_L1" "sigma_log_L2" "sigma_L2.1"  
```
