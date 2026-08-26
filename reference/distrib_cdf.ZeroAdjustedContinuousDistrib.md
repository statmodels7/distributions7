# Zero-Adjusted Continuous Cumulative Distribution Function

Computes \\F_Y(q) = (1-\pi)F_W(q;\theta) + \pi\\\mathbb{I}(q \ge 0)\\,
the parent's own distribution function scaled by \\1-\pi\\ with a JUMP
of size \\\pi\\ added at zero. That jump is the whole of what makes the
law mixed, and its height is exactly the probability
[`distrib_atoms.ZeroAdjustedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.ZeroAdjustedContinuousDistrib.md)
reports.

## Arguments

- distrib:

  A `ZeroAdjustedContinuousDistrib` object, from
  [`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md).

- q:

  A numeric vector of quantiles. The function is right continuous at
  zero, so `q = 0` includes the atom.

- theta:

  A named list with the parent's parameters followed by `za`.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, probabilities are \\P(Y
  \le q)\\; when `FALSE` they are \\P(Y \> q)\\.

- log.p:

  Logical of length 1. When `TRUE` the logarithm is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities, in \\\[0, 1\]\\.

## Notation

\\f_W\\ is the parent's density, \\\pi\\ the probability of the atom at
zero, \\f_Y\\ the mixed density and \\\ell\\ the log-density of one
observation.

## See also

[`distrib_pdf.ZeroAdjustedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.ZeroAdjustedContinuousDistrib.md)
for the density,
[`distrib_quantile.ZeroAdjustedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.ZeroAdjustedContinuousDistrib.md),
which inverts this across the jump, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- zero_adjusted(gaussian1_distrib())
theta <- list(mu = 1, sigma = 2, za = 0.3)

distrib_cdf(d, c(-1, 0, 2), theta)
#> [1] 0.1110587 0.5159763 0.7840237

# The jump at zero is the atom's probability.
distrib_cdf(d, 0, theta) - distrib_cdf(d, -1e-9, theta)
#> [1] 0.3

# Away from the atom it is the parent's, scaled and shifted.
all.equal(distrib_cdf(d, 2, theta), 0.7 * pnorm(2, 1, 2) + 0.3)
#> [1] TRUE
```
