# Expand Parameters to Common Length

Expands scalar parameters in a list to match the maximum length found
(or a specified length), ensuring all vectors are ready for element-wise
operations.

## Usage

``` r
expand_params(theta, n)
```

## Arguments

- theta:

  A named list of parameters.

- n:

  (Optional) The target length. If missing, defaults to
  `max(lengths(theta))`.

## Value

A list where all elements have length `n`.

## See also

[`deriv_names`](https://statmodels7.github.io/distributions7/reference/deriv_names.md),
[`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
[`transpose_params`](https://statmodels7.github.io/distributions7/reference/transpose_params.md),
[`check_params_dim`](https://statmodels7.github.io/distributions7/reference/check_params_dim.md),
[`check_theta_bounds`](https://statmodels7.github.io/distributions7/reference/check_theta_bounds.md),
[`param_smoothness`](https://statmodels7.github.io/distributions7/reference/param_smoothness.md),
[`generate_random_theta`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.md)

## Examples

``` r
expand_params(list(mu = 0, sigma = 1), n = 3)
#> $mu
#> [1] 0 0 0
#> 
#> $sigma
#> [1] 1 1 1
#> 
```
