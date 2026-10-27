# Expected Derivatives by Monte Carlo

Simulates `nsim` draws and averages over them the quantity of the
identity that is requested: the observed derivative for
`identity = "hessian"`, and the integrand of the score identity for
`identity = "score"`.

## Usage

``` r
expected_by_mc(distrib, y, theta, order, nsim, identity = "hessian")
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

- identity:

  `"hessian"` or `"score"`.

## Value

A named list of expected derivative component vectors; each has length
`length(y)`.

## Details

One simulation is run for each row of `theta` after recycling.
Observations whose parameters coincide are not merged, so a `theta` that
varies over the observations costs one simulation per observation, and a
scalar `theta` costs one simulation however long `y` is.

Estimates the same quantity as
[`expected_by_integrate()`](https://statmodels7.github.io/distributions7/reference/expected_by_integrate.md)
(Hessian identity) or
[`expected_by_bartlett()`](https://statmodels7.github.io/distributions7/reference/expected_by_bartlett.md)
(score identity), has an error that falls as \\1/\sqrt{n\_{sim}}\\, and
is stochastic unless the seed is fixed. Its cost is the cost of
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md),
so it is expensive for a distribution whose generator inverts a slow
quantile function. The Hessian identity at order 2 for a family with a
non-smooth parameter signals an error, as in
[`expected_by_integrate()`](https://statmodels7.github.io/distributions7/reference/expected_by_integrate.md).

## See also

[`expected_derivative_methods()`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
