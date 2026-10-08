# The Names of the Expected Information's Third Derivative

One key per pair \\(a,b)\\, spelled as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
spells it, and per unordered triple, spelled as
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
spells it at order three, joined `ab` first.

## Usage

``` r
d3expected_names(params)
```

## Arguments

- params:

  A character vector of parameter names, in the family's order.

## Value

A character vector,
`length(hess_names(params)) * length(deriv_names(params, 3))` long.

## See also

[`d3expected_key()`](https://statmodels7.github.io/distributions7/reference/d3expected_key.md),
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md)

## Examples

``` r
d3expected_names(c("mu", "sigma"))
#>  [1] "mu_mu_mu_mu_mu"                "mu_mu_mu_mu_sigma"            
#>  [3] "mu_mu_mu_sigma_sigma"          "mu_mu_sigma_sigma_sigma"      
#>  [5] "sigma_sigma_mu_mu_mu"          "sigma_sigma_mu_mu_sigma"      
#>  [7] "sigma_sigma_mu_sigma_sigma"    "sigma_sigma_sigma_sigma_sigma"
#>  [9] "mu_sigma_mu_mu_mu"             "mu_sigma_mu_mu_sigma"         
#> [11] "mu_sigma_mu_sigma_sigma"       "mu_sigma_sigma_sigma_sigma"   
```
