# Raw and Central Moments of a Distribution

Computes raw moments \\E\[Y^p\]\\ or central moments \\E\[(Y-\mu)^p\]\\
of a distribution numerically, via
[`expectation`](https://statmodels7.github.io/distributions7/reference/expectation.md)
(numerical integration for continuous distributions, series summation
for discrete ones).

## Usage

``` r
moment(distrib, theta, p = 1, central = FALSE, mu = NULL, ...)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- theta:

  A named list of parameters. Vectors are supported (vectorized
  computation).

- p:

  Numeric. The order of the moment. Can be a vector (recycled against
  `theta`).

- central:

  Logical. If `TRUE`, computes the central moment \\E\[(Y-\mu)^p\]\\,
  where \\\mu = E\[Y\]\\ (computed numerically unless `mu` is supplied).

- mu:

  Optional numeric. The centering value(s) used when `central = TRUE`.
  If `NULL`, the mean is computed numerically.

- ...:

  Additional arguments passed to
  [`expectation`](https://statmodels7.github.io/distributions7/reference/expectation.md).

## Value

A numeric vector of moments, with length equal to the maximum length
among `theta` components and `p`.

## See also

[`expectation`](https://statmodels7.github.io/distributions7/reference/expectation.md),
[`variance`](https://statmodels7.github.io/distributions7/reference/variance.md),
[`std_dev`](https://statmodels7.github.io/distributions7/reference/std_dev.md),
[`skewness`](https://statmodels7.github.io/distributions7/reference/skewness.md),
[`kurtosis`](https://statmodels7.github.io/distributions7/reference/kurtosis.md)

## Examples

``` r
if (FALSE) { # \dontrun{
d <- gaussian1_distrib()
moment(d, list(mu = 2, sigma = 3), p = 1)                 # 2
moment(d, list(mu = 2, sigma = 3), p = 2, central = TRUE) # 9
} # }
```
