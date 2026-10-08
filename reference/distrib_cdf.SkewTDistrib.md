# Skew t Distribution Function

Computes \\P(Y \le q)\\ by numerical integration of the density, taken
by the compiled rule of
[`compiled_cdf()`](https://statmodels7.github.io/distributions7/reference/compiled_cdf.md)
over the tail on the side of \\q\\ away from \\\mu\\.

## Arguments

- distrib:

  A `SkewTDistrib` object, from
  [`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu`, `sigma`, `alpha` and `nu`.

- lower.tail:

  Logical; if `FALSE`, \\P(Y \> q)\\ is returned.

- log.p:

  Logical; if `TRUE`, the logarithm is returned.

- ...:

  Unused.

## Value

A numeric vector of probabilities.

## Examples

``` r
d <- skewt_distrib()
th <- list(mu = 0.3, sigma = 1.2, alpha = 2, nu = 6)
distrib_cdf(d, c(-1, 0.5, 2), th)
#> [1] 0.004380413 0.219685163 0.795206582
c(method = distrib_cdf(d, 2, th),
  integral = integrate(function(v) distrib_pdf(d, v, th), -Inf, 2)$value)
#>    method  integral 
#> 0.7952066 0.7952066 
```
