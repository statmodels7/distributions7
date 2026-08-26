# Bernoulli Probability Mass Function

Computes the Bernoulli probability mass \$\$P(Y = y; \mu) =
\mu^{y}(1-\mu)^{1-y}, \qquad y \in \\0, 1\\,\$\$ by calling
[`stats::dbinom()`](https://rdrr.io/r/stats/Binomial.html) at
`size = 1`. The `pdf` in the generic's name is the density with respect
to counting measure, so what this returns is a probability.

## Arguments

- distrib:

  A `BernoulliDistrib` object, from
  [`bernoulli_distrib()`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md).

- y:

  A numeric vector of zeros and ones. Any other value gives 0 with a
  warning from
  [`stats::dbinom()`](https://rdrr.io/r/stats/Binomial.html).

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of the length of `y`. A value of length 1 is recycled.
  `mu` must lie in \\(0, 1)\\; a value outside \\\[0, 1\]\\ gives `NaN`
  with a warning.

- log:

  Logical of length 1. When `TRUE` the log-mass is returned, which stays
  finite for a probability so small that the mass underflows. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities in \\\[0, 1\]\\, of length
`max(length(y), length(mu))`, one value per observation.

## See also

[`distrib_pdf.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.BinomialDistrib.md)
for several trials,
[`distrib_gradient.BernoulliDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.BernoulliDistrib.md)
for the derivatives of the log-mass, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- bernoulli_distrib()

# The method is stats::dbinom at size 1.
all.equal(distrib_pdf(d, c(0, 1, 1), list(mu = 0.3)),
          dbinom(c(0, 1, 1), size = 1, prob = 0.3))
#> [1] TRUE

# The two masses are 1 - p and p, and sum to one.
distrib_pdf(d, c(0, 1), list(mu = 0.3))
#> [1] 0.7 0.3

# A probability may vary by observation, one value each.
distrib_pdf(d, c(0, 1, 1), list(mu = c(0.1, 0.5, 0.9)))
#> [1] 0.9 0.5 0.9
```
