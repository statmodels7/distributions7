# The Scale Matrix of a Multivariate Student t

Returns the scale matrix \\\Sigma\\, assembled from the matrix
parametrization's free values. This is the matrix in the density, and it
is not the covariance: the covariance is \\\nu\Sigma/(\nu-2)\\ and
exists only for \\\nu \> 2\\, while the scale matrix exists at every
admissible \\\nu\\.
[`variance.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.MvStudentTDistrib.md)
returns the covariance. The correlations read off either are the same, a
positive multiple of a matrix leaving them alone.

## Arguments

- distrib:

  An
  [MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object, from
  [`mvstudent_t_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t_distrib.md).

- theta:

  A named list of parameters, each component a single number. The
  location components and `nu` are ignored.

## Value

A \\p \times p\\ symmetric positive definite numeric matrix, with both
dimnames `v1`, ..., `vp`.

## See also

[`variance.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.MvStudentTDistrib.md)
for the covariance,
[`mv_sigma.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.MvGaussianDistrib.md),
where the two coincide,
[`mv_summary()`](https://statmodels7.github.io/distributions7/reference/mv_summary.md),
which reports the square roots of this matrix's diagonal as
`scale_sd_v1` and so on, and
[`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.md)
for the generic.

## Examples

``` r
d <- mvstudent_t_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)

mv_sigma(d, theta)
#>           v1        v2
#> v1 1.2214028 0.4420684
#> v2 0.4420684 0.8303200

# Positive definite whatever the free values are.
eigen(mv_sigma(d, theta), only.values = TRUE)$values
#> [1] 1.5092462 0.5424766

# It does not move with nu, where the covariance does.
vapply(c(3, 6, 60), function(nu) {
  t2 <- theta; t2$nu <- nu
  c(scale = mv_sigma(d, t2)[1, 1], covariance = variance(d, t2)[1, 1])
}, numeric(2))
#>                [,1]     [,2]     [,3]
#> scale      1.221403 1.221403 1.221403
#> covariance 3.664208 1.832104 1.263520
```
