# Transpose and Simplify Parameter List Structure

Transposes a list structure (swapping "columns" and "rows") and
simplifies the inner elements into atomic vectors.

Turns a list of `k` equal-length vectors into a list of `n` vectors of
length `k`, one per observation, keeping the names.

## Usage

``` r
transpose_params(theta)
```

## Arguments

- theta:

  A list to be transposed.

## Value

A `list` where each element has been transposed and simplified to an
atomic vector.

## See also

[`deriv_names`](https://statmodels7.github.io/distributions7/reference/deriv_names.md),
[`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
[`expand_params`](https://statmodels7.github.io/distributions7/reference/expand_params.md),
[`check_params_dim`](https://statmodels7.github.io/distributions7/reference/check_params_dim.md),
[`check_theta_bounds`](https://statmodels7.github.io/distributions7/reference/check_theta_bounds.md),
[`param_smoothness`](https://statmodels7.github.io/distributions7/reference/param_smoothness.md),
[`generate_random_theta`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.md)

## Examples

``` r
transpose_params(list(mu = c(0, 1), sigma = c(1, 2)))
#> [[1]]
#>    mu sigma 
#>     0     1 
#> 
#> [[2]]
#>    mu sigma 
#>     1     2 
#> 
```
