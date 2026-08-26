# Variance of the Generalized Gamma Distribution

Closed form from the first two raw moments, \\\operatorname{Var}(Y) =
m_2 - m_1^2\\ with \\m_k = a^k\\\Gamma\\(d+k)/p\\/\Gamma(d/p)\\. The
scale enters as a square and the two shapes through the gamma ratios.

## Arguments

- x:

  A `GenGamma1Distrib`, from
  [`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md).

- theta:

  A named list with components `a`, `d` and `p`, all positive, each a
  numeric vector of length 1 or `n`. Cancellation between the two raw
  moments costs digits where the coefficient of variation is small.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of length equal to the longest of the
three components.

## Notation

\\a \> 0\\ is the scale, \\d \> 0\\ and \\p \> 0\\ the two shapes, and
\\m_k\\ the \\k\\-th raw moment.

## See also

[`mean.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.GenGamma1Distrib.md),
[`skewness.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.GenGamma1Distrib.md),
[`gengamma_raw_moments()`](https://statmodels7.github.io/distributions7/reference/gengamma_raw_moments.md),
[`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md).

## Examples

``` r
d <- gengamma1_distrib()

# At p = 1 the family is a gamma of shape d, whose variance is a^2 d.
all.equal(variance(d, list(a = 2, d = 3, p = 1)), 12)
#> [1] TRUE

# At d = p it agrees with the Weibull.
all.equal(variance(d, list(a = 2, d = 3, p = 3)),
          variance(weibull1_distrib(), list(mu = 2, sigma = 3)))
#> [1] TRUE
```
