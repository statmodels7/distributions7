# Binomial Probability Mass Function

Computes the binomial probability mass \$\$P(Y = y; \mu) =
\binom{n}{y}\mu^{y}(1-\mu)^{n-y}, \qquad y = 0, 1, \dots, n,\$\$ by
calling [`stats::dbinom()`](https://rdrr.io/r/stats/Binomial.html) at
`size = distrib@size`. The number of trials comes from the object rather
than from `theta`, `size` being fixed data and not a parameter.

## Arguments

- distrib:

  A `BinomialDistrib` object, from
  [`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md).
  Its `size` property supplies the number of trials.

- y:

  A numeric vector of counts of successes, integers between 0 and
  `size`. Any other value gives 0 with a warning from
  [`stats::dbinom()`](https://rdrr.io/r/stats/Binomial.html).

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of the length of `y`. A value of length 1 is recycled.
  `mu` must lie in \\(0, 1)\\.

- log:

  Logical of length 1. When `TRUE` the log-mass is returned, which stays
  finite where the mass underflows. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities in \\\[0, 1\]\\, of length
`max(length(y), length(mu), length(distrib@size))`.

## See also

[`distrib_pdf.BernoulliDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.BernoulliDistrib.md)
for one trial,
[`distrib_cdf.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.BinomialDistrib.md)
for the distribution function, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- binomial_distrib(size = 10)

# The method is stats::dbinom at the object's own size.
all.equal(distrib_pdf(d, c(0, 4, 10), list(mu = 0.3)),
          dbinom(c(0, 4, 10), size = 10, prob = 0.3))
#> [1] TRUE

# A probability mass: it sums to one over 0:n.
sum(distrib_pdf(d, 0:10, list(mu = 0.3)))
#> [1] 1

# size may be one value per observation, which is what grouped binary data
# with unequal group sizes needs.
g <- binomial_distrib(size = c(5, 10, 20))
distrib_pdf(g, c(1, 4, 8), list(mu = 0.3))
#> [1] 0.3601500 0.2001209 0.1143967
```
