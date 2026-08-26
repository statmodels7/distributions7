# Generate Names for Higher-Order Derivative Components

Generates the names of the unique partial derivatives of a given `order`
with respect to a vector of parameters. Because mixed partial
derivatives are symmetric, only one representative per multi-index is
listed, using non-decreasing parameter order (e.g. `"mu_mu_sigma"` but
not `"mu_sigma_mu"`). For `order = 2` this coincides with the set of
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
(though possibly in a different order).

## Usage

``` r
deriv_names(params, order)
```

## Arguments

- params:

  A character vector of parameter names (e.g., `c("mu", "sigma")`).

- order:

  A positive integer, the derivative order (e.g. `3` or `4`).

## Value

A character vector of the \\\binom{p + \text{order} - 1}{\text{order}}\\
unique component names, where \\p\\ is the number of parameters.

## See also

[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
[`expand_params()`](https://statmodels7.github.io/distributions7/reference/expand_params.md),
[`transpose_params()`](https://statmodels7.github.io/distributions7/reference/transpose_params.md),
[`check_params_dim()`](https://statmodels7.github.io/distributions7/reference/check_params_dim.md),
[`check_theta_bounds()`](https://statmodels7.github.io/distributions7/reference/check_theta_bounds.md),
[`param_smoothness()`](https://statmodels7.github.io/distributions7/reference/param_smoothness.md),
[`generate_random_theta()`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.md)

## Examples

``` r
deriv_names(c("mu", "sigma"), 3)
#> [1] "mu_mu_mu"          "mu_mu_sigma"       "mu_sigma_sigma"   
#> [4] "sigma_sigma_sigma"
# "mu_mu_mu" "mu_mu_sigma" "mu_sigma_sigma" "sigma_sigma_sigma"
```
