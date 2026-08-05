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

## Examples

``` r
param_smoothness(gaussian1_distrib())
#>    mu sigma 
#>  TRUE  TRUE 

# the Laplace location is a kink, so it is not smooth
param_smoothness(laplace_distrib())
#>    mu     b 
#> FALSE  TRUE 
```
