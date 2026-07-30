# Map Parameters to the Link Scale

Applies each parameter's link, \\\eta_i = g_i(\theta_i)\\.

## Usage

``` r
fit_eta_from_theta(distrib, theta)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- theta:

  A named list of parameters on the natural scale.

## Value

A numeric vector of length `length(distrib@params)`.

## Details

Optimization is carried out on \\\eta\\, which is unconstrained, so this
is how a starting value expressed in natural parameters enters the
optimizer. Only the first element of each parameter is taken: a fit
estimates one \\\theta\\ for the whole sample.

## See also

[`fit_theta_from_eta`](https://statmodels7.github.io/distributions7/reference/fit_theta_from_eta.md),
the inverse.
