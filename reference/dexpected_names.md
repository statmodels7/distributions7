# The Names of the Expected Information's Derivative

One key per pair \\(a,b)\\ and differentiating parameter \\c\\, built by
joining
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
with the parameter differentiated in.

## Usage

``` r
dexpected_names(params)
```

## Arguments

- params:

  A character vector of parameter names, in the family's order.

## Value

A character vector, `length(hess_names(params)) * length(params)` long.

## Details

The keys are BUILT and never parsed, which is the package's rule
wherever a component name is a concatenation of parameter names: a
parameter whose own name contains an underscore makes the string
ambiguous to read back, and
[`dexpected_key()`](https://statmodels7.github.io/distributions7/reference/dexpected_key.md)
exists so that a consumer composes the same string this function
enumerates.

## See also

[`dexpected_key()`](https://statmodels7.github.io/distributions7/reference/dexpected_key.md),
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)

## Examples

``` r
dexpected_names(c("mu", "sigma"))
#> [1] "mu_mu_mu"          "mu_mu_sigma"       "sigma_sigma_mu"   
#> [4] "sigma_sigma_sigma" "mu_sigma_mu"       "mu_sigma_sigma"   
```
