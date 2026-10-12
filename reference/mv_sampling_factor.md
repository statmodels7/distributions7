# A Factor of the Covariance for Sampling

Returns a matrix \\F\\ with \\F F^\top = \Sigma\\, read from the matrix
parametrization's own factor instead of a Cholesky decomposition of the
assembled matrix. Where the parametrization carries \\\Sigma\\, \\F\\ is
its factor \\L\\; where it carries the precision \\\Omega = L L^\top\\,
\\F = L^{-\top}\\.

## Usage

``` r
mv_sampling_factor(s, eta, inverted)
```

## Arguments

- s:

  A parameters7 matrix parametrization of full rank.

- eta:

  Its free vector.

- inverted:

  Logical of length 1, `TRUE` when `s` carries the precision.

## Value

A \\p \times p\\ numeric matrix without dimnames.

## Details

A covariance that is nearly singular can round to a matrix that
[`chol()`](https://rdrr.io/r/base/chol.html) rejects although the
parametrization describes a valid one: at a log-Cholesky free value
\\\log L\_{22} = -20\\ the entry \\\Sigma\_{22} = 1 + e^{-40}\\ rounds
to 1. The factor of the parametrization is exact there.

## See also

[`distrib_rng.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.MvGaussianDistrib.md)
and
[`distrib_rng.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.MvStudentTDistrib.md),
the two callers, and
[`parameters7::param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.html)
for the factor.

## Examples

``` r
s <- parameters7::log_cholesky(2)
eta <- c(0, -20, 1)
f <- distributions7:::mv_sampling_factor(s, eta, FALSE)
f
#>      [,1]         [,2]
#> [1,]    1 0.000000e+00
#> [2,]    1 2.061154e-09
# chol() of the assembled matrix rejects it.
try(chol(parameters7::param_value(s, eta)))
#> Error in chol.default(parameters7::param_value(s, eta)) : 
#>   the leading minor of order 2 is not positive
```
