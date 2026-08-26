# Generate Names for Hessian Matrix Components

Generates the names of the unique second-order partial derivatives
(Hessian components) for a vector of parameter names: first the diagonal
elements (`"mu_mu"`, ...), then the upper-triangular off-diagonal
elements in row-major order (`"mu_sigma"`, ...).

## Usage

``` r
hess_names(params)
```

## Arguments

- params:

  A character vector of parameter names (e.g., `c("mu", "sigma")`).

## Value

A character vector of length \\n + n(n-1)/2\\.

## See also

[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md),
[`expand_params()`](https://statmodels7.github.io/distributions7/reference/expand_params.md),
[`transpose_params()`](https://statmodels7.github.io/distributions7/reference/transpose_params.md),
[`check_params_dim()`](https://statmodels7.github.io/distributions7/reference/check_params_dim.md),
[`check_theta_bounds()`](https://statmodels7.github.io/distributions7/reference/check_theta_bounds.md),
[`param_smoothness()`](https://statmodels7.github.io/distributions7/reference/param_smoothness.md),
[`generate_random_theta()`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.md)

## Examples

``` r
hess_names(c("mu", "sigma"))
#> [1] "mu_mu"       "sigma_sigma" "mu_sigma"   
# "mu_mu" "sigma_sigma" "mu_sigma"
```
