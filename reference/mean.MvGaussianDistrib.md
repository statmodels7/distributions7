# Mean of a Multivariate Gaussian

Returns \\\mathbb{E}\[Y\] = \mu\\, the mean vector. For this family the
expectation is a parameter and needs no integration: the density is
symmetric about \\\mu\\ and every coordinate has finite moments of every
order, so the first moment exists at every parameter value the
constructor admits.

## Arguments

- x:

  An
  [MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object, from
  [`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md).

- theta:

  A named list of parameters, each component a single number. Only the
  \\p\\ mean components are read.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length \\p\\, named `v1`, ..., `vp`.

## See also

[`variance.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.MvGaussianDistrib.md)
for the second moment,
[`mv_location.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_location.MvGaussianDistrib.md),
which returns the same vector as the density's center, and
[`base::mean()`](https://rdrr.io/r/base/mean.html) for the generic.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)

mean(d, theta)
#>   v1   v2 
#>  0.5 -0.3 

# Which is what a large sample average approaches.
set.seed(2)
round(colMeans(distrib_rng(d, 20000, theta)), 3)
#>     v1     v2 
#>  0.506 -0.299 
```
