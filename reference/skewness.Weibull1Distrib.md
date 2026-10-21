# Skewness of the Weibull Distribution

The skewness has the closed form \$\$\gamma_1 = \frac{g_3 - 3 g_1 g_2 +
2 g_1^3}{(g_2 - g_1^2)^{3/2}}, \qquad g_k = \Gamma(1 + k/\sigma).\$\$
The scale cancels, so the value depends only on the shape. It is zero at
a shape near 3.60235, so a Weibull is right-skewed below that value and
left-skewed above it, whereas the skewness of a gamma is always
positive. The method evaluates it through
[`weibull_moment_terms()`](https://statmodels7.github.io/distributions7/reference/weibull_moment_terms.md),
from the logarithms of the gamma factors and, for a shape of 5 or more,
from their series in \\1/\sigma\\, so that the cancellation among the
\\g_k\\ is never formed.

## Arguments

- x:

  A `Weibull1Distrib`, from
  [`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md).

- theta:

  A named list with components `mu` (the scale, positive) and `sigma`
  (the shape, positive); each is a numeric vector of length 1 or `n`.
  The third moment exists for every positive shape.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector, of length
`max(length(theta$mu), length(theta$sigma))`. Only the shape enters the
value, so a setting that varies only the scale repeats one number.

## Notation

\\\sigma \> 0\\ is the shape and \\g_k = \Gamma(1 + k/\sigma)\\.

## See also

[`kurtosis.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Weibull1Distrib.md),
also free of the scale;
[`variance.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Weibull1Distrib.md),
which is not;
[`weibull_moment_terms()`](https://statmodels7.github.io/distributions7/reference/weibull_moment_terms.md),
[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md).

## Examples

``` r
d <- weibull1_distrib()

# At shape 1 the family is exponential, whose skewness is 2.
all.equal(skewness(d, list(mu = 1, sigma = 1)), 2)
#> [1] TRUE

# The scale does not enter a standardized moment.
skewness(d, list(mu = c(0.1, 1, 100), sigma = 2))
#> [1] 0.6311107 0.6311107 0.6311107

# The sign changes at a shape of about 3.60235.
round(skewness(d, list(mu = 1, sigma = c(2, 3.60235, 10))), 6)
#> [1]  0.631111  0.000000 -0.637637
```
