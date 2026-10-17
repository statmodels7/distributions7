# Quantities Declared by the Matrix Parametrization

Returns the block that a parameters7 parametrization declares through
`param_readable()`, widened to the parameter vector of the distribution:
the Jacobian columns of the parametrization are placed at the positions
of its free values in `distrib@params`, and every other column is zero.
An AR(1) covariance, for instance, declares a scale and a correlation,
which a summary reports beside the standard deviations.

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
the shape of the return value of
[`mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.md):
`value`, `jacobian`, `transform` and `block`, with the Jacobian as wide
as `distrib@params`, and a logical attribute `precision`, `TRUE` when
the quantities describe a precision.

## Details

A parametrization that declares nothing gives `NULL`, and the summary is
then the result of
[`mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.md)
without it. The label ends in `(precision)` when the declared quantities
describe a precision, because the parametrization itself does not record
whether its matrix is a covariance or a precision. They describe a
precision when the distribution is inverted and the parametrization is
not an inverse one, or the other way round: an inverse parametrization
such as
[`parameters7::ar1_inv()`](https://statmodels7.github.io/parameters7/reference/ar1_inv.html)
declares, with `inverse = TRUE`, the quantities of the matrix that it
inverts, so under an inverted distribution they describe the covariance.

The positions of the free values are found by name: the family names
each free value with a prefix followed by the free name, after the
location for the gaussian and the Student t and first for the Dirichlet
and the multinomial. Where a name does not identify exactly one
parameter, the positions after the location are used. A multivariate
family written outside this package may have no `param` property, so the
function tests for it with
[`S7::prop_names()`](https://rconsortium.github.io/S7/reference/prop_names.html)
and returns `NULL` when it is absent.

## See also

[`parameters7::param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.html)
for what a parametrization declares,
[`mv_append_block()`](https://statmodels7.github.io/distributions7/reference/mv_append_block.md)
for how the result is joined on, and
[`mv_derived.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/mv_derived.multivariate_distrib.md)
and
[`mv_own_block()`](https://statmodels7.github.io/distributions7/reference/mv_own_block.md)
for the consumers.

## Examples

``` r
# A log-Cholesky covariance declares no block.
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0.5)
th <- distributions7:::align_theta(d, theta)
is.null(distributions7:::mv_param_block(d, th))
#> [1] TRUE

# An AR(1) covariance declares a scale and a correlation.
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

# The Jacobian is as wide as the parameter vector of the distribution, with
# the mean columns zero.
dim(pb$jacobian)
#> [1] 2 5
pb$jacobian[, 1:3]
#>       mu1 mu2 mu3
#> scale   0   0   0
#> rho     0   0   0

# The simplex of a Dirichlet declares its probabilities, whose free values
# come first among the parameters.
dd <- dirichlet_distrib(3)
thd <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8)
pbd <- distributions7:::mv_param_block(dd, distributions7:::align_theta(dd, thd))
pbd$value
#>        p1        p2        p3 
#> 0.4260125 0.2583897 0.3155978 
```
