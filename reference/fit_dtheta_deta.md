# Jacobian of the Inverse Link at the Estimate

The first derivative \\h_i'(\eta_i)\\ of each parameter's inverse link,
one entry per parameter.

## Usage

``` r
fit_dtheta_deta(distrib, eta)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- eta:

  A numeric vector of linear predictors, one per parameter.

## Value

A numeric vector of length `length(distrib@params)`.

## Details

This is the diagonal Jacobian the delta method needs to carry a standard
error from the link scale, where it is computed, to the parameter scale,
where it is reported.

## See also

[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
