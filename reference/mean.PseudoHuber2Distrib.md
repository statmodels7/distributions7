# Pseudo-Huber Moments, Standard-Deviation Parametrization

The mean is \\\mu\\, the variance \\\sigma^2\\ and the skewness zero.
The excess kurtosis depends on \\\nu\\ alone, \\3 K_3(\sqrt\nu)
K_1(\sqrt\nu) / K_2(\sqrt\nu)^2 - 3\\, which runs from 3 at \\\nu \to
0\\ to 0 at \\\nu \to \infty\\.

## Arguments

- x:

  A `PseudoHuber2Distrib` object.

- theta:

  A named list with components `mu`, `sigma` and `nu`.

- ...:

  Unused.

## Value

A numeric vector, one value per parameter row.

## See also

[`pseudohuber2_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber2_distrib.md).

## Examples

``` r
d <- pseudohuber2_distrib()
th <- list(mu = 0.4, sigma = 1.5, nu = 2)
c(mean(d, th), variance(d, th), skewness(d, th), kurtosis(d, th))
#> [1] 0.400000 2.250000 0.000000 1.534651
```
