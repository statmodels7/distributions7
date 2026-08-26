# Excess Kurtosis of the Generalized Gamma Distribution

Assembled from the first four raw moments, \\\gamma_2 = (m_4 - 4m_1m_3 +
6m_1^2m_2 - 3m_1^4)/(m_2 - m_1^2)^2 - 3\\ with \\m_k =
a^k\\\Gamma\\(d+k)/p\\/\Gamma(d/p)\\. The scale cancels, and the two
shapes move the skewness and the kurtosis with some freedom, where a
gamma ties them by \\\gamma_2 = 3\gamma_1^2/2\\.

## Arguments

- x:

  A `GenGamma1Distrib`, from
  [`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md).

- theta:

  A named list with components `a`, `d` and `p`, all positive, each a
  numeric vector of length 1 or `n`. Only `d` and `p` enter the value;
  the fourth-order combination cancels heavily at small dispersions.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of excess kurtoses, of length equal to the longest of
the three components.

## Notation

\\d \> 0\\ and \\p \> 0\\ are the two shapes and \\m_k\\ the \\k\\-th
raw moment.

## See also

[`skewness.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.GenGamma1Distrib.md),
from the same raw moments;
[`gengamma_raw_moments()`](https://statmodels7.github.io/distributions7/reference/gengamma_raw_moments.md),
[`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md).

## Examples

``` r
d <- gengamma1_distrib()

# At d = p it agrees with the Weibull.
all.equal(kurtosis(d, list(a = 2, d = 3, p = 3)),
          kurtosis(weibull1_distrib(), list(mu = 2, sigma = 3)))
#> [1] TRUE

# At p = 1 it agrees with the gamma, which ties it to the skewness.
th <- list(a = 2, d = 3, p = 1)
all.equal(kurtosis(d, th), 1.5 * skewness(d, th)^2)
#> [1] TRUE
```
