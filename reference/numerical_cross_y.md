# Numerical Mixed Response-Parameter Derivatives

Computes \\\partial^2 \ell / \partial y\\ \partial \theta_i\\ by one
central difference of
[`distrib_grad_y`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
in each parameter. Powers the default
[`distrib_cross_y`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
method for continuous distributions without a closed form.

## Usage

``` r
numerical_cross_y(distrib, y, theta, h_rel = .Machine$double.eps^(1/3))
```

## Arguments

- distrib:

  An object inheriting from class `"continuous_distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- h_rel:

  Numeric. Relative finite-difference step. Defaults to
  `.Machine$double.eps^(1/3)`.

## Value

A named list with one numeric vector per parameter.

## Details

The reference is the response gradient, not the log-density, so that a
distribution with an analytical `distrib_grad_y` pays for exactly one
finite-difference layer. When the response gradient is itself the
finite-difference fallback, the composition is the four-point mixed
stencil on the log-density – the two differences act on different
variables, so they commute into a single stencil rather than compounding
the way nested differences in the same variable do.

## See also

[`numerical_grad_y`](https://statmodels7.github.io/distributions7/reference/numerical_grad_y.md),
[`distrib_cross_y`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)

## Examples

``` r
numerical_cross_y(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#> $mu
#> [1] 1 1 1
#> 
#> $sigma
#> [1] -2  0  2
#> 
```
