# Generate Random Parameters

Generates sensible random parameters for a distribution based on its
mathematical domain.

## Usage

``` r
generate_random_theta(distrib, ...)
```

## Arguments

- distrib:

  A distribution object inheriting from the `distrib` class.

- ...:

  Additional arguments passed to the specific method.

## Value

A named list of parameter values, one per element of `distrib@params`,
each inside that parameter's bounds.

## See also

[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md),
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
[`expand_params()`](https://statmodels7.github.io/distributions7/reference/expand_params.md),
[`transpose_params()`](https://statmodels7.github.io/distributions7/reference/transpose_params.md),
[`check_params_dim()`](https://statmodels7.github.io/distributions7/reference/check_params_dim.md),
[`check_theta_bounds()`](https://statmodels7.github.io/distributions7/reference/check_theta_bounds.md),
[`param_smoothness()`](https://statmodels7.github.io/distributions7/reference/param_smoothness.md)

## Examples

``` r
set.seed(1)
generate_random_theta(gamma2_distrib())
#> $mu
#> [1] 1.666501
#> 
#> $sigma2
#> [1] 2.295531
#> 
```
