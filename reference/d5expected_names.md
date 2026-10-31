# The Names of the Expected Information's Fifth Derivative

Returns one key per pair \\(a,b)\\, spelled as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
spells it, and per unordered quintuple, spelled as
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
spells it at order five, joined `ab` first.

## Usage

``` r
d5expected_names(params)
```

## Arguments

- params:

  A character vector of parameter names, in the family's order.

## Value

A character vector,
`length(hess_names(params)) * length(deriv_names(params, 5))` long.

## See also

[`d5expected_key()`](https://statmodels7.github.io/distributions7/reference/d5expected_key.md),
[`d4expected_names()`](https://statmodels7.github.io/distributions7/reference/d4expected_names.md)

## Examples

``` r
d5expected_names(c("mu", "sigma"))
#>  [1] "mu_mu_mu_mu_mu_mu_mu"                     
#>  [2] "mu_mu_mu_mu_mu_mu_sigma"                  
#>  [3] "mu_mu_mu_mu_mu_sigma_sigma"               
#>  [4] "mu_mu_mu_mu_sigma_sigma_sigma"            
#>  [5] "mu_mu_mu_sigma_sigma_sigma_sigma"         
#>  [6] "mu_mu_sigma_sigma_sigma_sigma_sigma"      
#>  [7] "sigma_sigma_mu_mu_mu_mu_mu"               
#>  [8] "sigma_sigma_mu_mu_mu_mu_sigma"            
#>  [9] "sigma_sigma_mu_mu_mu_sigma_sigma"         
#> [10] "sigma_sigma_mu_mu_sigma_sigma_sigma"      
#> [11] "sigma_sigma_mu_sigma_sigma_sigma_sigma"   
#> [12] "sigma_sigma_sigma_sigma_sigma_sigma_sigma"
#> [13] "mu_sigma_mu_mu_mu_mu_mu"                  
#> [14] "mu_sigma_mu_mu_mu_mu_sigma"               
#> [15] "mu_sigma_mu_mu_mu_sigma_sigma"            
#> [16] "mu_sigma_mu_mu_sigma_sigma_sigma"         
#> [17] "mu_sigma_mu_sigma_sigma_sigma_sigma"      
#> [18] "mu_sigma_sigma_sigma_sigma_sigma_sigma"   
```
