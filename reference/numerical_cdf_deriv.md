# Numerical Derivatives of the Distribution Function

Central finite differences of
[`distrib_cdf`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
with respect to each parameter. These power the default
[`distrib_grad_cdf`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
and
[`distrib_hess_cdf`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md)
methods for continuous distributions that do not supply a closed form.

## Usage

``` r
numerical_cdf_deriv(
  distrib,
  q,
  theta,
  order = 1L,
  h_rel = .Machine$double.eps^(1/(order + 2)),
  which = NULL
)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters.

- order:

  Either 1 or 2.

- h_rel:

  Numeric. Relative finite-difference step.

- which:

  Character vector naming the components to differentiate, or `NULL`
  (default) for all of them: parameter names at first order and
  [`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
  components at second. Used by families that have a closed form for
  some components and not others, so that only the remaining ones cost
  cdf evaluations.

## Value

A named list of derivative components of \\F\\, not of its logarithm.

## See also

[`distrib_grad_cdf`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)

## Examples

``` r
numerical_cdf_deriv(gaussian1_distrib(), 1, list(mu = 0, sigma = 1), order = 1)
#> $mu
#> [1] -0.2419707
#> 
#> $sigma
#> [1] -0.2419707
#> 
```
