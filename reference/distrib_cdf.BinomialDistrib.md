# Binomial Cumulative Distribution Function

Computes the binomial distribution function \$\$F(q; \mu) = \sum\_{k =
0}^{\lfloor q \rfloor} \binom{n}{k}\mu^{k}(1-\mu)^{n-k}\$\$ by calling
[`stats::pbinom()`](https://rdrr.io/r/stats/Binomial.html) at
`size = distrib@size`, which evaluates it through the incomplete beta
function. The function is a step function, constant between integers,
and reaches 1 at `size`.

## Arguments

- distrib:

  A `BinomialDistrib` object, from
  [`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md).
  Its `size` property supplies the number of trials.

- q:

  A numeric vector of quantiles. A non-integer is floored, a negative
  value gives 0, and a value at or above `size` gives 1.

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of the length of `q`. `mu` must lie in \\(0, 1)\\.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, probabilities are \\P(Y
  \le q)\\; when `FALSE` they are \\P(Y \> q)\\, computed directly and
  so exact in the upper tail.

- log.p:

  Logical of length 1. When `TRUE` the logarithm of the probability is
  returned. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities in \\\[0, 1\]\\, of length
`max(length(q), length(mu), length(distrib@size))`.

## See also

[`distrib_quantile.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.BinomialDistrib.md)
for the generalized inverse,
[`distrib_pdf.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.BinomialDistrib.md)
for the mass, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- binomial_distrib(size = 10)
th <- list(mu = 0.3)

# The method is stats::pbinom at the object's own size.
all.equal(distrib_cdf(d, c(0, 4, 10), th),
          pbinom(c(0, 4, 10), size = 10, prob = 0.3))
#> [1] TRUE

# A step function reaching one at size.
distrib_cdf(d, c(3, 3.5, 3.9, 10, 20), th)
#> [1] 0.6496107 0.6496107 0.6496107 1.0000000 1.0000000

# The jump at an integer is exactly the mass there.
c(jump = distrib_cdf(d, 4, th) - distrib_cdf(d, 3, th),
  mass = distrib_pdf(d, 4, th))
#>      jump      mass 
#> 0.2001209 0.2001209 
```
