# Student t Cumulative Distribution Function

Computes the location-scale Student t distribution function \$\$F(q;
\mu, \sigma, \nu) = T\_\nu\\\left(\dfrac{q-\mu}{\sigma}\right)\$\$ with
\\T\_\nu\\ the standard Student t distribution function on \\\nu\\
degrees of freedom, by calling
[`stats::pt()`](https://rdrr.io/r/stats/TDist.html) at the standardized
value. Both tails are available exactly: `lower.tail = FALSE` evaluates
\\1 - F\\ without forming the difference, and `log.p = TRUE` returns a
logarithm that stays finite where the probability itself underflows.

## Arguments

- distrib:

  A `StudentT1Distrib` object, from
  [`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu`, `sigma` and `nu`, each a numeric
  vector of length 1 or of the length of `q`. A component of length 1 is
  recycled. `sigma` and `nu` must be strictly positive.

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
`max(length(q), length(mu), length(sigma), length(nu))`. With
`log.p = TRUE` the values are logarithms and are non-positive.

## See also

[`distrib_quantile.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.StudentT1Distrib.md)
for the inverse,
[`distrib_pdf.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.StudentT1Distrib.md)
for the density,
[`distrib_grad_cdf.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.StudentT1Distrib.md)
for the derivatives of this function in the parameters, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- student_t1_distrib()
th <- list(mu = 0.4, sigma = 1.2, nu = 5)

# The method is stats::pt at the standardized value.
all.equal(distrib_cdf(d, c(-2.5, 0.3, 1.8), th),
          pt((c(-2.5, 0.3, 1.8) - 0.4) / 1.2, df = 5))
#> [1] TRUE

# The law is symmetric about mu, so the two tails at equal distance match.
c(distrib_cdf(d, 0.4 - 2, th),
  distrib_cdf(d, 0.4 + 2, th, lower.tail = FALSE))
#> [1] 0.07822892 0.07822892

# The upper tail decays polynomially, so at forty scales out it is still
# representable where a Gaussian's has underflowed.
c(t = distrib_cdf(d, 48, list(mu = 0, sigma = 1.2, nu = 5),
                  lower.tail = FALSE),
  gaussian = pnorm(48, 0, 1.2, lower.tail = FALSE))
#>            t     gaussian 
#> 9.205981e-08 0.000000e+00 
```
