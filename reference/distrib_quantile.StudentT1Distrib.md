# Student t Quantile Function

Computes the location-scale Student t quantile function \$\$Q(p; \mu,
\sigma, \nu) = \mu + \sigma\\ T\_\nu^{-1}(p)\$\$ with \\T\_\nu^{-1}\\
the standard Student t quantile function on \\\nu\\ degrees of freedom,
by calling [`stats::qt()`](https://rdrr.io/r/stats/TDist.html). The
distribution function is strictly increasing on the whole line, so the
inverse is exact and unique and the root-finding fallback the base class
supplies is bypassed. `Q(0)` is `-Inf` and `Q(1)` is `Inf`.

## Arguments

- distrib:

  A `StudentT1Distrib` object, from
  [`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or of their
  logarithms when `log.p = TRUE`. A value outside the range gives `NaN`
  with a warning from
  [`stats::qt()`](https://rdrr.io/r/stats/TDist.html).

- theta:

  A named list with components `mu`, `sigma` and `nu`, each a numeric
  vector of length 1 or of the length of `p`. A component of length 1 is
  recycled. `sigma` and `nu` must be strictly positive.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, `p` is \\P(Y \le q)\\;
  when `FALSE` it is \\P(Y \> q)\\.

- log.p:

  Logical of length 1. When `TRUE`, `p` is read as a logarithm. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of quantiles on \\\[-\infty, \infty\]\\, of length
`max(length(p), length(mu), length(sigma), length(nu))`.

## See also

[`distrib_cdf.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.StudentT1Distrib.md)
for the function inverted here,
[`distrib_rng.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.StudentT1Distrib.md)
for draws, and
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- student_t1_distrib()
th <- list(mu = 0.4, sigma = 1.2, nu = 5)

# The quartiles, and the round trip back through the distribution function.
q <- distrib_quantile(d, c(0.25, 0.5, 0.75), th)
q
#> [1] -0.4720242  0.4000000  1.2720242
all.equal(distrib_cdf(d, q, th), c(0.25, 0.5, 0.75))
#> [1] TRUE

# The median is the location, the law being symmetric about it.
distrib_quantile(d, 0.5, th)
#> [1] 0.4

# Heavy tails put the extreme quantiles far further out than a Gaussian's.
rbind(t = distrib_quantile(d, c(0.001, 0.999), th),
      gaussian = qnorm(c(0.001, 0.999), 0.4, 1.2))
#>               [,1]     [,2]
#> t        -6.672115 7.472115
#> gaussian -3.308279 4.108279
```
