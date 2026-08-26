# Beta Cumulative Distribution Function in Mean and Precision

Computes the beta distribution function, the regularized incomplete beta
function \$\$F(q; \mu, \phi) = I_q(\alpha, \beta), \qquad \alpha =
\mu\phi, \quad \beta = (1-\mu)\phi,\$\$ by calling
[`stats::pbeta()`](https://rdrr.io/r/stats/Beta.html) at that pair of
shapes. Both tails are available exactly: `lower.tail = FALSE` evaluates
\\1 - F\\ without forming the difference, and `log.p = TRUE` returns a
logarithm that stays finite where the probability itself underflows to
zero.

## Arguments

- distrib:

  A `Beta1Distrib` object, from
  [`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md).

- q:

  A numeric vector of quantiles. A value at or below 0 gives a
  lower-tail probability of 0 and one at or above 1 gives 1.

- theta:

  A named list with components `mu` and `phi`, each a numeric vector of
  length 1 or of the length of `q`. A component of length 1 is recycled.
  `mu` must lie strictly in \\(0, 1)\\ and `phi` must be strictly
  positive.

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
`max(length(q), length(mu), length(phi))`. With `log.p = TRUE` the
values are logarithms and are non-positive.

## See also

[`distrib_quantile.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.Beta1Distrib.md)
for the inverse,
[`distrib_pdf.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Beta1Distrib.md)
for the density,
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
for the derivatives of this function in the parameters, which the beta
takes by finite difference because the derivative of an incomplete beta
in its shapes is hypergeometric, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- beta1_distrib()
th <- list(mu = 0.4, phi = 5)

# The method is stats::pbeta at the implied shapes.
all.equal(distrib_cdf(d, c(0.2, 0.5, 0.8), th),
          pbeta(c(0.2, 0.5, 0.8), 2, 3))
#> [1] TRUE

# The two tails sum to one.
distrib_cdf(d, 0.5, th) + distrib_cdf(d, 0.5, th, lower.tail = FALSE)
#> [1] 1

# At mu = 1/2 and phi = 2 the beta is the uniform, so F(q) = q.
distrib_cdf(d, c(0.2, 0.5, 0.8), list(mu = 0.5, phi = 2))
#> [1] 0.2 0.5 0.8

# Near the lower boundary the probability underflows and its log does not.
distrib_cdf(d, 1e-30, list(mu = 0.4, phi = 50))
#> [1] 0
distrib_cdf(d, 1e-30, list(mu = 0.4, phi = 50), log.p = TRUE)
#> [1] -1350.578
```
