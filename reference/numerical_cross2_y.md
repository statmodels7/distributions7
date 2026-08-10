# Numerical Mixed Second-Response Parameter Derivatives

Computes \\\partial^3 \ell / \partial y^2\\ \partial \theta_i\\ by one
central difference of
[`distrib_hess_y`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
in each parameter. Powers the default
[`distrib_cross2_y`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
method for continuous distributions without a closed form.

## Usage

``` r
numerical_cross2_y(
  distrib,
  y,
  theta,
  h_rel = .Machine$double.eps^(1/3),
  which = NULL
)
```

## Arguments

- distrib:

  A distribution object.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- h_rel:

  The relative step.

- which:

  Optional subset of parameters.

## Value

A named list with one numeric vector per parameter.

## Details

The reference is the response Hessian, not the log-density, so a
distribution with an analytical `distrib_hess_y` pays for exactly one
finite-difference layer. The two differences act on different variables
and therefore compose into a single mixed stencil rather than into the
nested differencing of one variable that the package forbids.

## Examples

``` r
numerical_cross2_y(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#> $mu
#> [1] 0 0 0
#> 
#> $sigma
#> [1] 2 2 2
#> 
```
