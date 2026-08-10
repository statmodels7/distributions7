# Expected Higher Derivatives of a Multivariate Family by Sampling

Averages the closed-form observed derivatives over a sample from the
family, which is the route the multivariate branch takes throughout: a
quadrature over the simplex has no counterpart to the one-dimensional
split at quantiles, and the Bartlett route would need the score's own
higher derivatives.

## Usage

``` r
mv_expected_higher(distrib, y, theta, order, nsim)
```

## Arguments

- distrib:

  A multivariate distribution object.

- y:

  The observed response, read only for its number of rows.

- theta:

  A named list of parameters.

- order:

  The derivative order, 3 or 4.

- nsim:

  The Monte Carlo sample size.

## Value

A named list of numeric vectors, each constant across observations.
