# Summed Score on the Link Scale

The gradient of the total log-likelihood with respect to \\\eta\\,
summed over observations.

## Usage

``` r
fit_score(distrib, y, theta, threads = 1L)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters on the natural scale.

## Value

A numeric vector of length `length(distrib@params)`.

## See also

[`fit_hess_matrix`](https://statmodels7.github.io/distributions7/reference/fit_hess_matrix.md)
