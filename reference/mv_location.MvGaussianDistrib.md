# Mean Vector of a Multivariate Gaussian

Returns the mean vector \\\mu\\, which for this family is the first
\\p\\ parameters read off `theta` in order. The method is
[`mv_leading_location()`](https://statmodels7.github.io/distributions7/reference/mv_leading_location.md),
shared with the Student t, whose location is also its leading
parameters. For the gaussian \\\mu\\ is both the center of the density
and the expectation of the law, which does not hold for every family.

## Arguments

- distrib:

  An
  [MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object, from
  [`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md).

- theta:

  A named list of parameters, one number each. Only the \\p\\ mean
  components are read.

## Value

A numeric vector of length \\p\\, named `v1`, ..., `vp` after the
coordinates of the response. The parameter names `mu1`, ..., `mup` do
not appear on it.

## See also

[`mean.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.MvGaussianDistrib.md),
which returns the same vector as the expectation of the law,
[`mv_sigma.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.MvGaussianDistrib.md)
for the matrix, and
[`mv_location()`](https://statmodels7.github.io/distributions7/reference/mv_location.md)
for the generic.

## Examples

``` r
d <- mvgaussian1_distrib(3)
theta <- as.list(stats::setNames(c(1, -2, 0.5, rep(0, 6)), d@params))

mv_location(d, theta)
#>   v1   v2   v3 
#>  1.0 -2.0  0.5 

# For a gaussian the location is also the expectation, so the two agree.
all.equal(unname(mv_location(d, theta)), unname(mean(d, theta)))
#> [1] TRUE
```
