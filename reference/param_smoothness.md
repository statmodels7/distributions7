# Per-Parameter Smoothness of the Log-Likelihood

Returns a named logical vector indicating, for each parameter of a
distribution, whether the log-likelihood is differentiable with respect
to it. This reads the `params_smooth` property, defaulting to all `TRUE`
when the property was left empty.

## Usage

``` r
param_smoothness(distrib)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

## Value

A named logical vector, one entry per parameter.

## See also

[`deriv_names`](https://statmodels7.github.io/distributions7/reference/deriv_names.md),
[`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
[`expand_params`](https://statmodels7.github.io/distributions7/reference/expand_params.md),
[`transpose_params`](https://statmodels7.github.io/distributions7/reference/transpose_params.md),
[`check_params_dim`](https://statmodels7.github.io/distributions7/reference/check_params_dim.md),
[`check_theta_bounds`](https://statmodels7.github.io/distributions7/reference/check_theta_bounds.md),
[`generate_random_theta`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.md)

## Examples

``` r
param_smoothness(gaussian1_distrib())
#>    mu sigma 
#>  TRUE  TRUE 

# the Laplace location is a kink, so it is not smooth
param_smoothness(laplace_distrib())
#>    mu sigma 
#> FALSE  TRUE 
```
