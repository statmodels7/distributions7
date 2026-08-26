# Bernoulli Cumulative Distribution Function

Computes the Bernoulli distribution function by calling
[`stats::pbinom()`](https://rdrr.io/r/stats/Binomial.html) at
`size = 1`. It takes three values only: \$\$F(q; \mu) = \begin{cases} 0,
& q \< 0,\\ 1 - \mu, & 0 \le q \< 1,\\ 1, & q \ge 1. \end{cases}\$\$

## Arguments

- distrib:

  A `BernoulliDistrib` object, from
  [`bernoulli_distrib()`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md).

- q:

  A numeric vector of quantiles. A non-integer is floored.

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of the length of `q`. `mu` must lie in \\(0, 1)\\.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, probabilities are \\P(Y
  \le q)\\; when `FALSE` they are \\P(Y \> q)\\.

- log.p:

  Logical of length 1. When `TRUE` the logarithm of the probability is
  returned. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities in \\\[0, 1\]\\, of length
`max(length(q), length(mu))`.

## See also

[`distrib_quantile.BernoulliDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.BernoulliDistrib.md)
for the generalized inverse,
[`distrib_pdf.BernoulliDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.BernoulliDistrib.md)
for the mass, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- bernoulli_distrib()
th <- list(mu = 0.3)

# Three values, whatever q is.
distrib_cdf(d, c(-1, 0, 0.5, 1, 2), th)
#> [1] 0.0 0.7 0.7 1.0 1.0

# The jump at 1 is exactly the mass there.
c(jump = distrib_cdf(d, 1, th) - distrib_cdf(d, 0, th),
  mass = distrib_pdf(d, 1, th))
#> jump mass 
#>  0.3  0.3 

# The method is stats::pbinom at size 1.
all.equal(distrib_cdf(d, c(0, 1), th), pbinom(c(0, 1), size = 1, prob = 0.3))
#> [1] TRUE
```
