# Mean of the Generalized Pareto Distribution

Closed form: \\E\[Y\] = \sigma/(1-\xi)\\ for \\\xi \< 1\\, and `Inf` at
or above one, where the first moment diverges. The shape controls which
moments exist at all: the moment of order \\k\\ is finite exactly for
\\\xi \< 1/k\\, so a fitted object can report a mean and no variance.

## Arguments

- x:

  A `GPDDistrib`, from
  [`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

- theta:

  A named list with components `sigma` (the scale, positive) and `xi`
  (the shape, any real value), each a numeric vector of length 1 or `n`.
  Settings with \\\xi \ge 1\\ give `Inf` in their own positions.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, of length
`max(length(theta$sigma), length(theta$xi))`, `Inf` wherever \\\xi \ge
1\\.

## Notation

\\\sigma \> 0\\ is the scale and \\\xi\\ the shape, of either sign. At
\\\xi = 0\\ the family is the exponential of mean \\\sigma\\, and at
\\\xi \< 0\\ the support is bounded above by \\\sigma/\|\xi\|\\.

## See also

[`variance.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.GPDDistrib.md),
whose threshold is \\1/2\\;
[`mean.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.ExponentialDistrib.md),
the \\\xi = 0\\ case;
[`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

## Examples

``` r
d <- gpd_distrib()

# Finite below shape one and infinite at or above it.
mean(d, list(sigma = 1, xi = c(0, 0.5, 1, 1.5)))
#> [1]   1   2 Inf Inf

# At shape zero the family is exponential of mean sigma.
all.equal(mean(d, list(sigma = 2, xi = 0)),
          mean(exponential_distrib(), list(mu = 2)))
#> [1] TRUE
```
