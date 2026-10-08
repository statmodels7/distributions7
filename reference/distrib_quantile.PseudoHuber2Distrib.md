# Pseudo-Huber Quantile Function, Standard-Deviation Parametrization

Inverts
[`distrib_cdf.PseudoHuber2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.PseudoHuber2Distrib.md)
by root finding on the lower half, in a bracket that starts ten standard
deviations below \\\mu\\ and is doubled until it contains the quantile;
an upper quantile is the reflection of the lower one about \\\mu\\.

## Arguments

- distrib:

  A `PseudoHuber2Distrib` object.

- p:

  A numeric vector of probabilities.

- theta:

  A named list with components `mu`, `sigma` and `nu`.

- lower.tail:

  Logical; if `FALSE`, `p` is an upper-tail probability.

- log.p:

  Logical; if `TRUE`, `p` is given on the log scale.

- ...:

  Unused.

## Value

A numeric vector of quantiles; `NaN` for a probability outside \\\[0,
1\]\\.

## See also

[`distrib_cdf.PseudoHuber2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.PseudoHuber2Distrib.md).

## Examples

``` r
d <- pseudohuber2_distrib()
th <- list(mu = 0.4, sigma = 1.5, nu = 2)
q <- distrib_quantile(d, c(0.1, 0.5, 0.9), th)
distrib_cdf(d, q, th)
#> [1] 0.1 0.5 0.9
```
