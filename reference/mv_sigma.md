# The Matrix That a Parameter List Describes

Returns the matrix that describes the spread of a multivariate
distribution, as a \\p \times p\\ numeric matrix: the covariance for the
gaussian, the Dirichlet and the multinomial, and the scale matrix for
the Student t. The Student t's covariance is \\\nu\Sigma/(\nu-2)\\ for
\\\nu \> 2\\ and is infinite for \\\nu \le 2\\, while its scale matrix
exists at every \\\nu\\;
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
is the generic for the second moment.

## Usage

``` r
mv_sigma(distrib, theta)
```

## Arguments

- distrib:

  An object inheriting from
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md).

- theta:

  A named list or vector of parameters, each component a single number.
  Aligned by the generic before dispatch.

## Value

A \\p \times p\\ symmetric numeric matrix: positive definite, with both
dimnames `v1`, ..., `vp`, for the gaussian and the Student t; of rank
\\p - 1\\ and without dimnames for the Dirichlet and the multinomial,
whose coordinates sum to a constant.

## Details

Where the parametrization carries the precision, as for
[`mvgaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md)
and
[`mvstudent_t2_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md),
the matrix is inverted here, so the result is the covariance or the
scale matrix either way.

No method is registered on the class, so a family without a method of
its own fails in dispatch.

## Notation

\\\Sigma\\ is the matrix that the parametrization carries, \\\nu\\ the
degrees of freedom of a Student t and \\p\\ the dimension.

## See also

[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
for the moment,
[`mv_location()`](https://statmodels7.github.io/distributions7/reference/mv_location.md)
for the location,
[`mv_summary()`](https://statmodels7.github.io/distributions7/reference/mv_summary.md)
for the standard deviations and correlations of a fit, and
[`mv_sigma.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.MvGaussianDistrib.md)
and
[`mv_sigma.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.MvStudentTDistrib.md)
for two methods.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 1, mu2 = -1, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0.5)
mv_sigma(d, theta)
#>     v1   v2
#> v1 1.0 0.50
#> v2 0.5 1.25

# For a gaussian the matrix is the variance; for a Student t it is not.
t2 <- mvstudent_t1_distrib(2)
th <- c(theta, list(nu = 6))
all.equal(mv_sigma(d, theta), variance(d, theta))
#> [1] TRUE
all.equal(variance(t2, th), (6 / 4) * mv_sigma(t2, th))
#> [1] TRUE

# At two degrees of freedom the covariance is infinite and the scale matrix
# is finite.
th2 <- th; th2$nu <- 2
c(scale = mv_sigma(t2, th2)[1, 1], covariance = variance(t2, th2)[1, 1])
#>      scale covariance 
#>          1        Inf 
```
