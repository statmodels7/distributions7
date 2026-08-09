# Quantile Function

Evaluates the quantile function for a given distribution.

## Usage

``` r
distrib_quantile(distrib, p, theta, ...)
```

## Arguments

- distrib:

  A distribution object inheriting from the `distrib` class.

- p:

  A numeric vector of probabilities.

- theta:

  A named list (or named numeric vector) of distribution parameters. If
  unnamed, parameters are taken in the order of `distrib@params`.

- ...:

  Additional arguments passed to the specific method (e.g.,
  `lower.tail`, `log.p`).

## Value

A numeric vector of quantiles.

## Details

The generalized inverse of the distribution function,

\$\$Q(p; \theta) = \inf\\y : F(y; \theta) \ge p\\,\$\$

which for a continuous strictly increasing \\F\\ is the ordinary inverse
and for a discrete family the smallest support point whose cumulative
mass reaches \\p\\. Without a method the value comes from root-finding
on
[`distrib_cdf`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
in the continuous case and from an exact table lookup in the discrete
one.

## See also

[`distrib_pdf`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md),
[`distrib_cdf`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md),
[`distrib_rng`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md),
[`distrib_atoms`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md)

## Examples

``` r
distrib_quantile(gaussian1_distrib(), c(0.025, 0.5, 0.975), list(mu = 0, sigma = 1))
#> [1] -1.959964  0.000000  1.959964
distrib_quantile(poisson_distrib(), c(0.1, 0.9), list(mu = 2))
#> [1] 0 4
```
