# Folded Distribution Function

Computes \\P(\|Y\| \le q) = F(q) - F(-q)\\ from the parent's own
distribution function, exactly and with no quadrature. It is the
difference of two calls on the parent, clamped to \\\[0, 1\]\\ against
rounding and set to zero below the support.

## Arguments

- distrib:

  A `FoldedDistrib` object, from
  [`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md).

- q:

  A numeric vector of quantiles. Values below zero give probability `0`.

- theta:

  A named list of the parent's parameters.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, probabilities are
  \\P(\|Y\| \le q)\\; when `FALSE` they are \\P(\|Y\| \> q)\\.

- log.p:

  Logical of length 1. When `TRUE` the logarithm is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities, in \\\[0, 1\]\\.

## Details

`lower.tail = FALSE` and `log.p = TRUE` are formed from the computed
probability, as \\1 - p\\ and \\\log p\\, not by a separate route
through the parent. Far into the upper tail that difference cancels, so
a caller who needs the survival function of a folded distribution to
many digits should expect the loss the subtraction implies.

## Notation

\\F\\ is the parent's distribution function and \\Y\\ the parent's
variable.

## See also

[`distrib_pdf.FoldedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.FoldedDistrib.md)
for the density,
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md),
which the class inherits and which inverts this exactly, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- folded(gaussian1_distrib())
theta <- list(mu = 0.5, sigma = 1.2)

q <- c(0, 1, 3)
distrib_cdf(d, q, theta)
#> [1] 0.0000000 0.5558891 0.9796206

# Which is the difference of the parent's own distribution function.
all.equal(distrib_cdf(d, q, theta),
          pnorm(q, 0.5, 1.2) - pnorm(-q, 0.5, 1.2))
#> [1] TRUE

# Both tails and the logarithm.
distrib_cdf(d, 1, theta, lower.tail = FALSE)
#> [1] 0.4441109
distrib_cdf(d, 1, theta, log.p = TRUE)
#> [1] -0.5871865

# The inherited quantile inverts it, so the round trip closes.
p <- c(0.1, 0.5, 0.9)
max(abs(distrib_cdf(d, distrib_quantile(d, p, theta), theta) - p))
#> [1] 5.364487e-12
```
