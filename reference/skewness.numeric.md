# Sample Skewness

The sample skewness of a numeric vector, \\m_3 / m_2^{3/2}\\ with \\m_k
= n^{-1}\sum_i (y_i - \bar y)^k\\. Both central moments divide by \\n\\,
so this is the population denominator and the estimator is biased
towards zero in small samples.
[`variance.numeric()`](https://statmodels7.github.io/distributions7/reference/variance.numeric.md)
uses \\n - 1\\ instead, each following the commonest convention for its
own quantity.

## Arguments

- x:

  A numeric vector. A constant vector gives `NaN`, the standardizing
  denominator being zero.

- na.rm:

  Should missing values be dropped before the moments are taken? A
  single logical, `FALSE` by default, which propagates `NA`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A single number.

## See also

[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md)
for the generic and the distribution version,
[`kurtosis.numeric()`](https://statmodels7.github.io/distributions7/reference/kurtosis.numeric.md)
for the fourth order,
[`variance.numeric()`](https://statmodels7.github.io/distributions7/reference/variance.numeric.md),
whose denominator is \\n - 1\\.

## Examples

``` r
set.seed(1)
y <- rgamma(500, shape = 2)

# The n denominator, written out.
m <- mean(y)
all.equal(skewness(y),
          mean((y - m)^3) / mean((y - m)^2)^1.5)
#> [1] TRUE

# It estimates the population value 2 / sqrt(shape).
c(sample = skewness(y), theory = 2 / sqrt(2))
#>   sample   theory 
#> 1.385948 1.414214 

# Missing values propagate unless dropped.
c(skewness(c(y, NA)), skewness(c(y, NA), na.rm = TRUE))
#> [1]       NA 1.385948
```
