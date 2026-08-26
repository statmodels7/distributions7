# Elastic-Net Quantile Function

Inverts
[`distrib_cdf.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.EnetDistrib.md)
in closed form, each half of the density being a truncated Gaussian.
Below the median the quantile is \\\mu + \\x +
\Phi^{-1}(2p\\\Phi(-x))\\/\sqrt c\\ and above it the reflection. Nothing
is inverted by root finding, so
[`distrib_rng.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.EnetDistrib.md)
can use the inverse transform.

The `qnorm` call is made on the **log** scale, its argument being
\\\log(2p) + \log\Phi(-x)\\, because \\\Phi(-x)\\ underflows to zero
past \\x = 38\\.

## Arguments

- distrib:

  An `EnetDistrib` object, from
  [`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md).

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or their logarithms
  when `log.p = TRUE`.

- theta:

  A named list with components `mu`, `lambda` and `alpha`, each a
  numeric vector of length 1 or of the length of `p`.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, `p` is \\P(Y \le q)\\;
  when `FALSE` it is the survival probability.

- log.p:

  Logical of length 1. When `TRUE`, `p` is given as a logarithm and is
  exponentiated first. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of quantiles, of the length of the recycled inputs.

## See also

[`distrib_cdf.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.EnetDistrib.md),
which it inverts,
[`distrib_rng.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.EnetDistrib.md),
which draws from it, and
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- enet_distrib()
th <- list(mu = 0, lambda = 2, alpha = 0.5)

# The round trip is exact, both functions being closed form.
q <- c(-2, -0.5, 0.5, 2)
all.equal(distrib_quantile(d, distrib_cdf(d, q, th), th), q)
#> [1] TRUE

# Symmetric about the location.
p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
round(distrib_quantile(d, p, th), 10)
#> [1] -0.8559401 -0.4096087  0.0000000  0.4096087  0.8559401

# It still answers where the Mills argument has passed 38.
distrib_quantile(d, c(0.25, 0.75),
                 list(mu = 0, lambda = 20, alpha = 0.995))
#> [1] -0.03481968  0.03481968
```
