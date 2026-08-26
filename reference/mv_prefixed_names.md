# Prefix a Parametrization's Free Names with the Matrix They Describe

Returns the matrix parametrization's free names with `"sigma_"` or
`"omega_"` in front, so that a printed parameter table says which matrix
the coordinate belongs to. A free value's name says how the matrix is
BUILT, not which matrix it is, so the same parametrization on the two
sides of a model would otherwise give two genuinely different models the
same parameter names.

## Usage

``` r
mv_prefixed_names(free_names, inverted = FALSE)
```

## Arguments

- free_names:

  The parametrization's own free names, a character vector.

- inverted:

  Logical of length 1. `TRUE` for a precision, giving `"omega_"`;
  `FALSE`, the default, for a covariance or a scale matrix, giving
  `"sigma_"`.

## Value

A character vector as long as `free_names`.

## Details

The prefix is applied by the DISTRIBUTION. The parametrization does not
know which side of a model it has been handed to, so it cannot apply
one.

## See also

[`mvgaussian_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian_distrib.md),
whose two forms this distinguishes, and
[`parameters7::log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.html)
for a source of free names.

## Examples

``` r
distributions7:::mv_prefixed_names(c("log_L1", "log_L2", "L2.1"))
#> [1] "sigma_log_L1" "sigma_log_L2" "sigma_L2.1"  
distributions7:::mv_prefixed_names(c("log_L1", "L2.1"), inverted = TRUE)
#> [1] "omega_log_L1" "omega_L2.1"  

# Which is what separates the two parametrizations of one law.
mvgaussian_distrib(2)@params
#> [1] "mu1"          "mu2"          "sigma_log_L1" "sigma_log_L2" "sigma_L2.1"  
mvgaussian_distrib(2, omega = parameters7::log_cholesky(2))@params
#> [1] "mu1"          "mu2"          "omega_log_L1" "omega_log_L2" "omega_L2.1"  
```
