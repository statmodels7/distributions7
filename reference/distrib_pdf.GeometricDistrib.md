# Geometric Probability Mass Function

Computes the geometric probability mass \$\$P(Y = y; \mu) =
\dfrac{1}{1+\mu}\left(\dfrac{\mu}{1+\mu}\right)^{y}, \qquad y = 0, 1, 2,
\dots\$\$ by calling
[`stats::dgeom()`](https://rdrr.io/r/stats/Geometric.html) at
`prob = 1/(1+mu)`, which
[`geom_prob()`](https://statmodels7.github.io/distributions7/reference/geom_prob.md)
supplies. The mass falls geometrically, by the constant factor
\\\mu/(1+\mu)\\ at every step, so its maximum is always at zero however
large the mean is.

## Arguments

- distrib:

  A `GeometricDistrib` object, from
  [`geometric_distrib()`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md).

- y:

  A numeric vector of counts. A non-integer or negative value gives 0
  with a warning from
  [`stats::dgeom()`](https://rdrr.io/r/stats/Geometric.html).

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of the length of `y`. A value of length 1 is recycled.
  `mu` must be strictly positive.

- log:

  Logical of length 1. When `TRUE` the log-mass is returned, which stays
  finite for a count far above the mean. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities in \\\[0, 1\]\\, of length
`max(length(y), length(mu))`, one value per observation.

## See also

[`geom_prob()`](https://statmodels7.github.io/distributions7/reference/geom_prob.md)
for the conversion,
[`distrib_cdf.GeometricDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.GeometricDistrib.md)
for the distribution function, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- geometric_distrib()
y <- c(0, 2, 7)

# The method is stats::dgeom at prob = 1/(1+mu).
all.equal(distrib_pdf(d, y, list(mu = 3)), dgeom(y, prob = 1 / (1 + 3)))
#> [1] TRUE

# A probability mass: it sums to one over the support.
sum(distrib_pdf(d, 0:400, list(mu = 3)))
#> [1] 1

# The ratio of consecutive masses is the constant mu/(1+mu), so the mode is
# always at zero.
m <- distrib_pdf(d, 0:5, list(mu = 3))
c(ratios = unique(round(m[-1] / m[-length(m)], 12)), mu_over_1_plus_mu = 3 / 4)
#>            ratios mu_over_1_plus_mu 
#>              0.75              0.75 
```
