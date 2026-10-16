# Elastic-Net Quantile Function

Inverts
[`distrib_cdf.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.EnetDistrib.md)
in closed form, because each half of the density is a truncated
Gaussian. Below the median the quantile is \\\mu + \\x +
\Phi^{-1}(2p\\\Phi(-x))\\/\sqrt c\\ and above it the reflection. The
inversion does not use root finding, so
[`distrib_rng.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.EnetDistrib.md)
can use the inverse transform.

The probability is first turned into the logarithms of both tails,
without forming \\1 - p\\ or exponentiating a logarithm, and each half
is inverted from the logarithm of its own tail. The inversion calls
`qnorm` on the log scale, where its argument is \\\log P + \log 2 +
\log\Phi(-x)\\, because \\\Phi(-x)\\ underflows to zero beyond \\x
\approx 38.5\\. A probability outside \\\[0, 1\]\\ gives `NaN`.

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

  Logical of length 1. When `TRUE`, `p` is given as a logarithm.
  Defaults to `FALSE`.

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

# The round trip is exact because both functions are closed form.
q <- c(-2, -0.5, 0.5, 2)
all.equal(distrib_quantile(d, distrib_cdf(d, q, th), th), q)
#> [1] TRUE

# The quantiles are symmetric about the location.
p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
round(distrib_quantile(d, p, th), 10)
#> [1] -0.8559401 -0.4096087  0.0000000  0.4096087  0.8559401

# The quantile function still returns values where the Mills argument
# exceeds 38.
distrib_quantile(d, c(0.25, 0.75),
                 list(mu = 0, lambda = 20, alpha = 0.995))
#> [1] -0.03481968  0.03481968
```
