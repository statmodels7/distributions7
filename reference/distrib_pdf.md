# Probability Density Function

Evaluates the probability density function (PDF) or probability mass
function (PMF).

## Usage

``` r
distrib_pdf(distrib, y, theta, ...)
```

## Arguments

- distrib:

  A distribution object inheriting from the `distrib` class.

- y:

  A numeric vector of observations.

- theta:

  A named list (or named numeric vector) of distribution parameters. If
  unnamed, parameters are taken in the order of `distrib@params`.

- ...:

  Additional arguments passed to the specific method (e.g., `log`).

## Value

A numeric vector of density values, one per observation.

## Details

For a continuous family \\f(y; \theta)\\ is the density with respect to
the Lebesgue measure and for a discrete one the mass \\f(y; \theta) =
P(Y = y)\\; `log = TRUE` returns \\\log f(y; \theta)\\, which is the
quantity every derivative generic differentiates. This is the only
method a distribution has to supply: every other quantity has a
numerical fallback derived from it.

## See also

[`distrib_cdf`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md),
[`distrib_quantile`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md),
[`distrib_rng`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md),
[`distrib_atoms`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md)

## Examples

``` r
distrib_pdf(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#> [1] 0.2419707 0.3989423 0.2419707
distrib_pdf(poisson_distrib(), 0:3, list(mu = 2), log = TRUE)
#> [1] -2.000000 -1.306853 -1.306853 -1.712318
```
