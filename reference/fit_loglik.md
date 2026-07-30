# Total Log-Likelihood

The sum of the log-density over the observations, the objective the fit
maximises.

## Usage

``` r
fit_loglik(distrib, y, theta)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters on the natural scale.

## Value

A single number.

## See also

[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
