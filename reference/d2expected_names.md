# The Names of the Expected Information's Second Derivative

One key per pair \\(a,b)\\ and pair \\(c,d)\\, each pair spelled as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
spells it, joined `ab` first.

## Usage

``` r
d2expected_names(params)
```

## Arguments

- params:

  A character vector of parameter names, in the family's order.

## Value

A character vector, `length(hess_names(params))^2` long.

## See also

[`d2expected_key()`](https://statmodels7.github.io/distributions7/reference/d2expected_key.md),
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)

## Examples

``` r
d2expected_names(c("mu", "sigma"))
#> [1] "mu_mu_mu_mu"             "mu_mu_sigma_sigma"      
#> [3] "mu_mu_mu_sigma"          "sigma_sigma_mu_mu"      
#> [5] "sigma_sigma_sigma_sigma" "sigma_sigma_mu_sigma"   
#> [7] "mu_sigma_mu_mu"          "mu_sigma_sigma_sigma"   
#> [9] "mu_sigma_mu_sigma"      
```
