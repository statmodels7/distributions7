# The Quantities the Matrix Parametrization Is About

Returns the block a parameters7 family declares through
`param_readable()`, widened to the distribution's own parameter vector:
the parametrization's Jacobian columns are placed in the stretch of
`distrib@params` its free values occupy, and every other column is zero.
An AR(1) covariance, for instance, is about a scale and a correlation,
and those are what a reader wants beside the standard deviations.

## Usage

``` r
mv_param_block(distrib, theta)
```

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- theta:

  A named list of parameters, already aligned.

## Value

`NULL` where the parametrization declares nothing, or a named list in
the shape of
[`mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.md)'s
return: `value`, `jacobian`, `transform` and `block`, with the Jacobian
as wide as `distrib@params`.

## Details

A family that declares nothing returns `NULL`, and the summary is then
whatever
[`mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.md)
produced without it. The label says `(precision)` when the distribution
is inverted, because the parametrization does not know which side of the
model it was handed to.

A multivariate family written OUTSIDE this package may have no `param`
property at all, so the property is asked for with
[`S7::prop_names()`](https://rconsortium.github.io/S7/reference/prop_names.html)
rather than assumed.

## See also

[`parameters7::param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.html)
for what a parametrization declares,
[`mv_append_block()`](https://statmodels7.github.io/distributions7/reference/mv_append_block.md)
for how the result is joined on, and
[`mv_derived.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_derived.MvGaussianDistrib.md)
for the consumer.

## Examples

``` r
# A log-Cholesky covariance declares nothing: its free values are
# coordinates and nothing more.
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0.5)
th <- distributions7:::align_theta(d, theta)
is.null(distributions7:::mv_param_block(d, th))
#> [1] TRUE

# An AR(1) covariance is about a scale and a correlation, and says so.
a <- mvgaussian1_distrib(3, parameters7::ar1(3))
th <- as.list(stats::setNames(c(0, 0, 0, 0.1, 0.3), a@params))
pb <- distributions7:::mv_param_block(a, distributions7:::align_theta(a, th))
pb$value
#>     scale       rho 
#> 1.1051709 0.2913126 
pb$transform
#>   scale     rho 
#>   "log" "atanh" 
pb$block
#>                      scale                        rho 
#> "Autoregressive structure" "Autoregressive structure" 

# The Jacobian is as wide as the distribution's parameter vector, with the
# mean columns zero.
dim(pb$jacobian)
#> [1] 2 5
pb$jacobian[, 1:3]
#>       mu1 mu2 mu3
#> scale   0   0   0
#> rho     0   0   0
```
