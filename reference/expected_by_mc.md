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

A named list of expected derivative component vectors, each of length
`length(y)`.

## Details

One simulation is run per *distinct* parameter configuration rather than
per observation: the expectation depends on \\\theta\\ alone, so
observations sharing a \\\theta\\ share an answer, and a model with a
scalar \\\theta\\ costs one simulation however long `y` is.

Estimates the same quantity as
[`expected_by_integrate`](https://statmodels7.github.io/distributions7/reference/expected_by_integrate.md),
with an error falling as \\1/\sqrt{n\_{sim}}\\, and is stochastic unless
the seed is fixed. Its cost is the cost of sampling, so it is the wrong
choice for a distribution whose RNG falls back to inverse transform on a
numerical quantile.

## See also

[`expected_derivative_methods`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
