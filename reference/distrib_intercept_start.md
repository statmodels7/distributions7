# Starting Values for a Regression's Intercepts, Read Off the Data

Returns, on the **parameter** scale, a starting value for each parameter
whose intercept a regression should not take from the intercept-only
maximum likelihood fit. A modelling layer carries each value onto the
parameter's own link and writes it to that equation's intercept. The
base method returns an empty list: for most families the intercept-only
fit is the right start.

## Usage

``` r
distrib_intercept_start(distrib, y, ...)
```

## Arguments

- distrib:

  An object inheriting from `distrib`.

- y:

  The response, a numeric vector.

- ...:

  Passed to methods. No shipped method reads it.

## Value

A named list of single numbers on the parameter scale, possibly empty,
each strictly inside its parameter's bounds and named after an entry of
`distrib@params`.

## Details

The case this exists for is
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md).
Without covariates the parent's overdispersion can absorb the excess
zeros, so the intercept-only maximum puts the mixing weight at the edge
of its domain – measured, \\\pi = 2.2 \times 10^{-308}\\ on a negative
binomial with a true \\\pi = 0.25\\ – which is right for that model and
a trap as a start for a model with covariates: the link is flat there,
the score on the unconstrained scale vanishes, and the fit stays at the
edge reporting convergence while an interior maximum 7.6 log-likelihood
units higher exists. See
[`distrib_intercept_start.ZeroInflatedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_intercept_start.ZeroInflatedDistrib.md).

## See also

[`distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.md),
the starting values
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
uses.

## Examples

``` r
distrib_intercept_start(poisson_distrib(), c(0, 1, 3))
#> list()
distrib_intercept_start(zero_inflated(poisson_distrib()), c(0, 0, 1, 3))
#> $zi
#> [1] 0.5
#> 
```
