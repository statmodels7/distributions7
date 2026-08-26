# Bernoulli Log-CDF Gradient

Closed form, the binomial identity at \\n = 1\\. The distribution
function takes two values, \\F(0) = 1 - p\\ and \\F(1) = 1\\, so the
derivative is \\-1\\ at \\k = 0\\ and exactly zero at \\k = 1\\, the
upper value being 1 whatever \\p\\ is.

## Arguments

- distrib:

  A `BernoulliDistrib` object, from
  [`bernoulli_distrib()`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md).

- q:

  A numeric vector of quantiles. Non-integer values are floored; values
  below zero give a derivative of zero.

- theta:

  A named list with one component, `mu` (the success probability,
  strictly between 0 and 1), a numeric vector of length 1 or `n`.

- lower.tail:

  Is the lower tail wanted? A single logical, `TRUE` by default.

- log:

  Are derivatives of the log probability wanted? A single logical,
  `TRUE` by default. At \\k = 1\\ the lower-tail probability is 1 and
  the log-scale derivative is 0; the upper-tail probability is 0 there,
  and the log-scale derivative is `NaN`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with one numeric vector, `mu`, the length of `q` recycled
against `theta`.

## Notation

\\p \in (0,1)\\ is the success probability and \\k = \lfloor q
\rfloor\\.

## See also

[`distrib_grad_cdf.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.BinomialDistrib.md),
of which this is the \\n = 1\\ case;
[`bernoulli_distrib()`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md).

## Examples

``` r
d <- bernoulli_distrib()

# -1 at zero, and exactly 0 at one.
distrib_grad_cdf(d, c(0, 1), list(mu = 0.3), log = FALSE)$mu
#> [1] -1  0

# On the log scale at k = 0 that is -1 / (1 - p).
all.equal(distrib_grad_cdf(d, 0, list(mu = 0.3))$mu, -1 / 0.7)
#> [1] TRUE
```
