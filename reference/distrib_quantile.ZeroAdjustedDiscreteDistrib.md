# Zero-Adjusted Discrete Quantile Function

Inverts the hurdle distribution function. The quantile is `0` for \\p
\le \pi\\, the whole mass at zero, and otherwise the PARENT's quantile
at \\u\\1 - f(0)\\ + f(0)\\ with \\u = (p - \pi)/(1 - \pi)\\. The parent
is inverted in its upper tail, at \\(1 - p)\\1 - f(0)\\/(1 - \pi)\\, at
the rows of the observations, so that a constant that varies by
observation (a binomial's size) is read at its own row.

## Arguments

- distrib:

  A `ZeroAdjustedDiscreteDistrib` object, from
  [`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md).

- p:

  A numeric vector of probabilities, clamped to \\\[0, 1\]\\ after the
  `log.p` and `lower.tail` transformations are applied.

- theta:

  A named list with the parent's parameters followed by `za`.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, `p` is \\P(Y \le q)\\;
  when `FALSE` it is \\P(Y \> q)\\.

- log.p:

  Logical of length 1. When `TRUE`, `p` is given as a logarithm.
  Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of quantiles, of the recycled length of `p` and
`theta`. Every value is either `0` or a positive support point of the
parent.

## Details

The result is a quantile of a discrete distribution: the distribution
function evaluated at it is at least `p`, and usually larger. Every
probability at or below \\\pi\\ maps to zero, so the first step has
probability \\\pi\\, and `p = 1` gives the upper end of the support.

## Notation

\\f\\ is the parent's mass function, \\F\\ its distribution function,
\\\pi\\ the probability of a zero, \\s\\ the parent's score and \\\ell\\
the log-mass of one observation.

## See also

[`distrib_cdf.ZeroAdjustedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.ZeroAdjustedDiscreteDistrib.md),
which this inverts, and
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- zero_adjusted(poisson_distrib())
theta <- list(mu = 3, za = 0.4)

distrib_quantile(d, c(0.2, 0.5, 0.95), theta)
#> [1] 0 2 6

# Everything at or below the mass at zero maps to zero.
distrib_quantile(d, c(0.1, 0.4), theta)
#> [1] 0 0

# The round trip overshoots, as on any discrete support.
p <- c(0.5, 0.95)
rbind(asked = p,
      reached = distrib_cdf(d, distrib_quantile(d, p, theta), theta))
#>              [,1]      [,2]
#> asked   0.5000000 0.9500000
#> reached 0.6357806 0.9788415
```
