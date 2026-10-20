# Prefix Free Names with the Matrix That a Parametrization Describes

Returns the free names of a matrix parametrization with `"sigma_"` or
`"omega_"` in front, so that a printed parameter table names the matrix
to which each coordinate belongs. The free name of a value records how
the matrix is built and not which matrix it is, so without the prefix
the same parametrization of a covariance and of a precision would give
two different models the same parameter names.

## Usage

``` r
mv_prefixed_names(free_names, inverted = FALSE)
```

## Arguments

- free_names:

  The parametrization's own free names, a character vector.

- inverted:

  Logical of length 1. It is `TRUE` for a precision, which gives
  `"omega_"`, and `FALSE` (the default) for a covariance or a scale
  matrix, which gives `"sigma_"`.

## Value

A character vector as long as `free_names`.

## Details

The distribution applies the prefix, because the parametrization does
not record which side of a model it describes.

## See also

[`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md),
whose two forms this distinguishes, and
[`parameters7::log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.html)
for a source of free names.

## Examples

``` r
distributions7:::mv_prefixed_names(c("log_L1", "log_L2", "L2.1"))
#> [1] "sigma_log_L1" "sigma_log_L2" "sigma_L2.1"  
distributions7:::mv_prefixed_names(c("log_L1", "L2.1"), inverted = TRUE)
#> [1] "omega_log_L1" "omega_L2.1"  

# The prefixes separate the two parametrizations of one law.
mvgaussian1_distrib(2)@params
#> [1] "mu1"          "mu2"          "sigma_log_L1" "sigma_log_L2" "sigma_L2.1"  
mvgaussian2_distrib(2, parameters7::log_cholesky(2))@params
#> [1] "mu1"          "mu2"          "omega_log_L1" "omega_log_L2" "omega_L2.1"  
```
