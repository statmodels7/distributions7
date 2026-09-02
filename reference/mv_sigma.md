# The Matrix a Parameter List Describes

Returns the matrix a multivariate distribution's parametrization
carries, as a \\p \times p\\ numeric matrix: the COVARIANCE for a
gaussian, and the SCALE MATRIX for a Student t. It is not in general the
second moment. The Student t's covariance is \\\nu\Sigma/(\nu-2)\\ and
does not exist below two degrees of freedom, while its scale matrix
exists at every \\\nu\\;
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
is the generic that answers about the moment, and keeping the two apart
is what allows a heavy-tailed family to be described at all.

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

A \\p \times p\\ symmetric positive definite numeric matrix, with both
dimnames `v1`, ..., `vp`.

## Details

Where the parametrization carries the PRECISION, as
[`mvgaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md)
does, the matrix is inverted here, so the result is the covariance
either way.

The base-class method signals an error: not every multivariate family
has a matrix parameter, and a family that does registers its own method.

## Notation

\\\Sigma\\ is the matrix the parametrization carries, \\\nu\\ the
degrees of freedom of a Student t and \\p\\ the dimension.

## See also

[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
for the moment,
[`mv_location()`](https://statmodels7.github.io/distributions7/reference/mv_location.md)
for the location,
[`mv_summary()`](https://statmodels7.github.io/distributions7/reference/mv_summary.md)
for the standard deviations and correlations a reader wants, and
[`mv_sigma.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.MvGaussianDistrib.md)
and
[`mv_sigma.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.MvStudentTDistrib.md)
for the two methods.

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

# At two degrees of freedom the covariance is gone and the scale matrix
# stands.
th2 <- th; th2$nu <- 2
c(scale = mv_sigma(t2, th2)[1, 1], covariance = variance(t2, th2)[1, 1])
#>      scale covariance 
#>          1        Inf 
```
