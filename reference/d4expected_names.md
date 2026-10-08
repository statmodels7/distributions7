# The Names of the Expected Information's Fourth Derivative

One key per pair \\(a,b)\\, spelled as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
spells it, and per unordered quadruple, spelled as
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
spells it at order four, joined `ab` first.

## Usage

``` r
d4expected_names(params)
```

## Arguments

- params:

  A character vector of parameter names, in the family's order.

## Value

A character vector,
`length(hess_names(params)) * length(deriv_names(params, 4))` long.

## See also

[`d4expected_key()`](https://statmodels7.github.io/distributions7/reference/d4expected_key.md),
[`d3expected_names()`](https://statmodels7.github.io/distributions7/reference/d3expected_names.md)

## Examples

``` r
d4expected_names(c("mu", "sigma"))
#>  [1] "mu_mu_mu_mu_mu_mu"                   "mu_mu_mu_mu_mu_sigma"               
#>  [3] "mu_mu_mu_mu_sigma_sigma"             "mu_mu_mu_sigma_sigma_sigma"         
#>  [5] "mu_mu_sigma_sigma_sigma_sigma"       "sigma_sigma_mu_mu_mu_mu"            
#>  [7] "sigma_sigma_mu_mu_mu_sigma"          "sigma_sigma_mu_mu_sigma_sigma"      
#>  [9] "sigma_sigma_mu_sigma_sigma_sigma"    "sigma_sigma_sigma_sigma_sigma_sigma"
#> [11] "mu_sigma_mu_mu_mu_mu"                "mu_sigma_mu_mu_mu_sigma"            
#> [13] "mu_sigma_mu_mu_sigma_sigma"          "mu_sigma_mu_sigma_sigma_sigma"      
#> [15] "mu_sigma_sigma_sigma_sigma_sigma"   
```
