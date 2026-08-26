# Variance of a Multivariate Gaussian

Returns \\\operatorname{Var}(Y) = \Sigma\\, the covariance matrix the
matrix parametrization carries, inverted first where that
parametrization carries the precision. For this family the scale matrix
of the density is the covariance of the law, so the value agrees with
[`mv_sigma.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.MvGaussianDistrib.md);
the two part company in the heavy-tailed sibling, where the scale matrix
exists at every \\\nu\\ and the covariance does not.

The return is a matrix, not the numeric vector a univariate family
gives, the second central moment of a vector being a matrix.

## Arguments

- x:

  An
  [MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object, from
  [`mvgaussian_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian_distrib.md).

- theta:

  A named list of parameters, each component a single number. The \\p\\
  mean components are ignored.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A \\p \times p\\ symmetric positive definite numeric matrix, with both
dimnames `v1`, ..., `vp`.

## See also

[`mean.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.MvGaussianDistrib.md)
for the first moment,
[`mv_sigma.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.MvGaussianDistrib.md)
for the same matrix read as the density's scale,
[`variance.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.MvStudentTDistrib.md)
for the family where the two differ, and
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
for the generic.

## Examples

``` r
d <- mvgaussian_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)

variance(d, theta)
#>           v1        v2
#> v1 1.2214028 0.4420684
#> v2 0.4420684 0.8303200

# The scale matrix and the covariance are one matrix for this family.
identical(variance(d, theta), mv_sigma(d, theta))
#> [1] TRUE

# Which a large sample covariance approaches.
set.seed(2)
round(var(distrib_rng(d, 20000, theta)), 3)
#>       v1    v2
#> v1 1.236 0.451
#> v2 0.451 0.839
```
