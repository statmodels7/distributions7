# Geometric Cumulative Distribution Function

Computes the geometric distribution function \$\$F(q; \mu) = 1 -
\left(\dfrac{\mu}{1+\mu}\right)^{\lfloor q \rfloor + 1}\$\$ by calling
[`stats::pgeom()`](https://rdrr.io/r/stats/Geometric.html) at
`prob = 1/(1+mu)`. The survival function is exactly
\\(\mu/(1+\mu))^{q+1}\\, a geometric decay, and `lower.tail = FALSE`
returns it without forming the difference.

## Arguments

- distrib:

  A `GeometricDistrib` object, from
  [`geometric_distrib()`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md).

- q:

  A numeric vector of quantiles. A non-integer is floored, and a
  negative value gives 0.

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of the length of `q`. `mu` must be strictly positive.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, probabilities are \\P(Y
  \le q)\\; when `FALSE` they are the survival function \\P(Y \> q)\\,
  exact far into the tail.

- log.p:

  Logical of length 1. When `TRUE` the logarithm of the probability is
  returned. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities in \\\[0, 1\]\\, of length
`max(length(q), length(mu))`.

## See also

[`distrib_quantile.GeometricDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.GeometricDistrib.md)
for the generalized inverse,
[`distrib_pdf.GeometricDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.GeometricDistrib.md)
for the mass, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- geometric_distrib()
th <- list(mu = 3)

# The method is stats::pgeom at prob = 1/(1+mu).
all.equal(distrib_cdf(d, c(0, 2, 7), th), pgeom(c(0, 2, 7), prob = 1 / 4))
#> [1] TRUE

# The survival function is a geometric decay, written out.
all.equal(distrib_cdf(d, c(0, 2, 7), th, lower.tail = FALSE),
          (3 / 4)^(c(0, 2, 7) + 1))
#> [1] TRUE

# Memoryless: the chance of at least one more failure does not depend on
# how many have already been seen.
vapply(c(0, 2, 10, 50), function(s)
  distrib_cdf(d, s + 1, th, lower.tail = FALSE) /
    distrib_cdf(d, s, th, lower.tail = FALSE), numeric(1))
#> [1] 0.75 0.75 0.75 0.75
```
