# Expected Derivatives by Monte Carlo

Simulates `nsim` draws and averages the observed derivative over them.

## Usage

``` r
expected_by_mc(distrib, y, theta, order, nsim)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations; only its length is used.

- theta:

  A named list of parameters.

- order:

  The derivative order, 2 to 4.

- nsim:

  The number of draws per parameter configuration.

## Value

A named list of expected derivative component vectors; each has length
`length(y)`.

## Details

One simulation is run for each row of `theta` after recycling.
Observations whose parameters coincide are not merged, so a `theta` that
varies over the observations costs one simulation per observation, and a
scalar `theta` costs one simulation however long `y` is.

Estimates the same quantity as
[`expected_by_integrate()`](https://statmodels7.github.io/distributions7/reference/expected_by_integrate.md),
has an error that falls as \\1/\sqrt{n\_{sim}}\\, and is stochastic
unless the seed is fixed. Its cost is the cost of
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md),
so it is expensive for a distribution whose generator inverts a slow
quantile function.

## See also

[`expected_derivative_methods()`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
