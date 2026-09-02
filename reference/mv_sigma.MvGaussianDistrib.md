# Covariance of a Multivariate Gaussian

Returns the covariance \\\Sigma\\, assembled from the matrix
parametrization's free values. Where the parametrization carries the
precision the matrix is inverted first, so the result is the covariance
either way. For this family the matrix is both the scale matrix of the
density and the variance of the law, so
[`variance.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.MvGaussianDistrib.md)
returns the same matrix; the two part company in the heavy-tailed
sibling.

## Arguments

- distrib:

  An
  [MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object, from
  [`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md).

- theta:

  A named list of parameters, one number each. The \\p\\ mean components
  are ignored.

## Value

A \\p \times p\\ symmetric positive definite numeric matrix, with both
dimnames `v1`, ..., `vp`.

## See also

[`variance.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.MvGaussianDistrib.md)
for the same matrix read as a moment,
[`mv_location.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_location.MvGaussianDistrib.md)
for the mean,
[`mv_summary()`](https://statmodels7.github.io/distributions7/reference/mv_summary.md)
for the standard deviations and correlations a reader wants, and
[`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.md)
for the generic.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
round(mv_sigma(d, theta), 4)
#>        v1     v2
#> v1 1.2214 0.4421
#> v2 0.4421 0.8303

# For a gaussian the scale matrix is the variance.
all.equal(mv_sigma(d, theta), variance(d, theta))
#> [1] TRUE

# A precision parametrization reports the covariance too, inverted here.
o <- mvgaussian2_distrib(2, parameters7::log_cholesky(2))
th_o <- list(mu1 = 0, mu2 = 0, omega_log_L1 = 0.1,
             omega_log_L2 = -0.2, omega_L2.1 = 0.4)
Om <- parameters7::param_value(o@param, unlist(th_o)[3:5])
all.equal(mv_sigma(o, th_o), solve(Om), check.attributes = FALSE)
#> [1] TRUE
```
