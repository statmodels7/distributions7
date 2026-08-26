# Beta Cumulative Distribution Function in the Shapes

Computes the beta distribution function, the regularized incomplete beta
function \\F(q) = I_q(\alpha, \beta)\\, by calling
[`stats::pbeta()`](https://rdrr.io/r/stats/Beta.html) at
`shape1 = alpha` and `shape2 = beta`. Both tails are available exactly:
`lower.tail = FALSE` evaluates \\1 - F\\ without forming the difference,
and `log.p = TRUE` returns a logarithm that stays finite where the
probability itself underflows to zero.

## Arguments

- distrib:

  A `Beta2Distrib` object, from
  [`beta2_distrib()`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md).

- q:

  A numeric vector of quantiles. A value at or below 0 gives a
  lower-tail probability of 0 and one at or above 1 gives 1.

- theta:

  A named list with components `alpha` and `beta`, each a numeric vector
  of length 1 or of the length of `q`. A component of length 1 is
  recycled. Both must be strictly positive.

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
`max(length(q), length(alpha), length(beta))`. With `log.p = TRUE` the
values are logarithms and are non-positive.

## See also

[`distrib_quantile.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.Beta2Distrib.md)
for the inverse,
[`distrib_pdf.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Beta2Distrib.md)
for the density,
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
for the derivatives of this function in the parameters, which the beta
takes by finite difference because the derivative of an incomplete beta
in its shapes is hypergeometric, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- beta2_distrib()
th <- list(alpha = 2, beta = 5)

# The method is stats::pbeta at these two shapes.
all.equal(distrib_cdf(d, c(0.1, 0.3, 0.7), th),
          pbeta(c(0.1, 0.3, 0.7), 2, 5))
#> [1] TRUE

# The two tails sum to one.
distrib_cdf(d, 0.3, th) + distrib_cdf(d, 0.3, th, lower.tail = FALSE)
#> [1] 1

# Both shapes 1 is the uniform, so F(q) = q.
distrib_cdf(d, c(0.1, 0.3, 0.7), list(alpha = 1, beta = 1))
#> [1] 0.1 0.3 0.7

# Near the lower boundary the probability underflows and its log does not.
distrib_cdf(d, 1e-30, list(alpha = 20, beta = 5))
#> [1] 0
distrib_cdf(d, 1e-30, list(alpha = 20, beta = 5), log.p = TRUE)
#> [1] -1372.28
```
