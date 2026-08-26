# Poisson Cumulative Distribution Function

Computes the Poisson distribution function \$\$F(q; \mu) = P(Y \le q) =
\sum\_{k = 0}^{\lfloor q \rfloor} \dfrac{e^{-\mu}\mu^{k}}{k!}\$\$ by
calling [`stats::ppois()`](https://rdrr.io/r/stats/Poisson.html), which
evaluates it through the incomplete gamma function. The function is a
step function, constant between integers, so `F(2)` and `F(2.9)` are the
same number.

## Arguments

- distrib:

  A `PoissonDistrib` object, from
  [`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md).

- q:

  A numeric vector of quantiles. A non-integer is floored, and a
  negative value gives 0.

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of the length of `q`. `mu` must be strictly positive.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, probabilities are \\P(Y
  \le q)\\; when `FALSE` they are \\P(Y \> q)\\, computed directly and
  so exact far into the upper tail.

- log.p:

  Logical of length 1. When `TRUE` the logarithm of the probability is
  returned. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities in \\\[0, 1\]\\, of length
`max(length(q), length(mu))`.

## See also

[`distrib_quantile.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.PoissonDistrib.md)
for the generalized inverse,
[`distrib_pdf.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.PoissonDistrib.md)
for the mass, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- poisson_distrib()
th <- list(mu = 3)

# The method is stats::ppois.
all.equal(distrib_cdf(d, c(0, 2, 7), th), ppois(c(0, 2, 7), lambda = 3))
#> [1] TRUE

# A step function: constant between integers.
distrib_cdf(d, c(2, 2.5, 2.9), th)
#> [1] 0.4231901 0.4231901 0.4231901

# The jump at an integer is exactly the mass there.
c(jump = distrib_cdf(d, 2, th) - distrib_cdf(d, 1, th),
  mass = distrib_pdf(d, 2, th))
#>      jump      mass 
#> 0.2240418 0.2240418 

# The upper tail is exact where its complement would have rounded to one.
distrib_cdf(d, 60, th, lower.tail = FALSE)
#> [1] 1.310782e-56
```
