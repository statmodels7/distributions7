# Skewness of the Generalized Gamma Distribution

Assembled from the first three raw moments, \\\gamma_1 = (m_3 -
3m_1m_2 + 2m_1^3)/(m_2 - m_1^2)^{3/2}\\ with \\m_k =
a^k\\\Gamma\\(d+k)/p\\/\Gamma(d/p)\\. The scale cancels, so the value
depends on the two shapes alone; unlike a gamma's it can be negative,
which is part of what the second shape parameter buys.

## Arguments

- x:

  A `GenGamma1Distrib`, from
  [`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md).

- theta:

  A named list with components `a`, `d` and `p`, all positive, each a
  numeric vector of length 1 or `n`. Only `d` and `p` enter the value.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector, of length equal to the longest of the three
components.

## Notation

\\d \> 0\\ and \\p \> 0\\ are the two shapes and \\m_k\\ the \\k\\-th
raw moment. The scale does not enter a standardized moment.

## See also

[`kurtosis.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.GenGamma1Distrib.md),
from the same raw moments;
[`gengamma_raw_moments()`](https://statmodels7.github.io/distributions7/reference/gengamma_raw_moments.md),
[`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md).

## Examples

``` r
d <- gengamma1_distrib()

# At d = p it agrees with the Weibull, which changes sign with its shape.
all.equal(skewness(d, list(a = 2, d = 3, p = 3)),
          skewness(weibull1_distrib(), list(mu = 2, sigma = 3)))
#> [1] TRUE

# The scale does not enter a standardized moment.
skewness(d, list(a = c(0.1, 1, 100), d = 3, p = 1.5))
#> [1] 0.7375295 0.7375295 0.7375295
```
