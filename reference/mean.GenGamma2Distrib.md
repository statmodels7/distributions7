# Generalized Gamma Moments, Mean

With \\g_j = \Gamma((d+j)/p)\\\Gamma(d/p)^{j-1}/\Gamma((d+1)/p)^j\\, the
ratio \\E\[Y^j\]/E\[Y\]^j\\, the mean is \\m\\, the variance \\m^2
(g_2 - 1)\\, the skewness \\(g_3 - 3 g_2 + 2)/(g_2 - 1)^{3/2}\\ and the
excess kurtosis \\(g_4 - 4 g_3 + 6 g_2 - 3)/(g_2 - 1)^2 - 3\\.

## Arguments

- x:

  A `GenGamma2Distrib` object.

- theta:

  A named list with components `mean`, `d` and `p`.

- ...:

  Unused.

## Value

A numeric vector, one value per parameter row.

## Details

These combinations of the ratios are not formed: towards the lognormal
they cancel to several orders. The central moments of \\Y/m\\ and the
fourth cumulant come from the series of the cumulant generating function
of \\\log Y\\, in double-double arithmetic, as
[`variance.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.GenGamma1Distrib.md)
describes; the skewness is \\\mu_3/\mu_2^{3/2}\\ and the excess kurtosis
\\(\mu_4 - 3\mu_2^2)/\mu_2^2\\.

## See also

[`gengamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma2_distrib.md).

## Examples

``` r
d <- gengamma2_distrib()
th <- list(mean = 5, d = 3, p = 1.5)
c(mean(d, th), variance(d, th), skewness(d, th), kurtosis(d, th))
#> [1] 5.0000000 5.6809682 0.7375295 0.6362894
```
