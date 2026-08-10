# Numerical Hyperparameter Hessians of the Response Derivatives

Computes the second-order mixed components by one central difference of
an analytic first-order component in each parameter.

## Usage

``` r
numerical_theta2_y(distrib, y, theta, inner, h_rel = .Machine$double.eps^(1/3))
```

## Arguments

- distrib:

  A distribution object.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- inner:

  A function of `theta` returning the first-order components, one per
  parameter.

- h_rel:

  The relative step.

## Value

A named list keyed by `hess_names(distrib@params)`.

## Details

The reference is
[`distrib_cross_y`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
or
[`distrib_cross2_y`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md),
so a distribution with a closed form for those pays for exactly one
difference. A mixed pair is differenced both ways and averaged: the two
agree in exact arithmetic and not quite in floating point, the steps
differing, and a second derivative of a scalar has to come out
symmetric.

## Examples

``` r
d <- gaussian1_distrib()
numerical_theta2_y(d, c(-1, 0, 1), list(mu = 0, sigma = 1),
                   function(th) distrib_cross_y(d, c(-1, 0, 1), th))
#> $mu_mu
#> [1] 0 0 0
#> 
#> $sigma_sigma
#> [1]  6  0 -6
#> 
#> $mu_sigma
#> [1] -2 -2 -2
#> 
```
