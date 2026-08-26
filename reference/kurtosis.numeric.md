# Sample Excess Kurtosis

The sample excess kurtosis of a numeric vector, \\m_4 / m_2^{2} - 3\\
with \\m_k = n^{-1}\sum_i (y_i - \bar y)^k\\. Both central moments
divide by \\n\\, so this is the population denominator, the convention
[`skewness.numeric()`](https://statmodels7.github.io/distributions7/reference/skewness.numeric.md)
also uses. The estimator is biased downwards in small samples, and heavy
tails make it slow to settle: it depends on the fourth moment, so a
handful of extreme points dominate it.

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

[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md)
for the generic and the distribution version,
[`skewness.numeric()`](https://statmodels7.github.io/distributions7/reference/skewness.numeric.md)
for the third order,
[`variance.numeric()`](https://statmodels7.github.io/distributions7/reference/variance.numeric.md),
whose denominator is \\n - 1\\.

## Examples

``` r
set.seed(1)
y <- rt(1000, df = 10)

# The n denominator, written out.
m <- mean(y)
all.equal(kurtosis(y), mean((y - m)^4) / mean((y - m)^2)^2 - 3)
#> [1] TRUE

# It estimates 6 / (nu - 4), and needs a lot of data to get there.
c(sample = kurtosis(y), theory = 6 / (10 - 4))
#>   sample   theory 
#> 1.130882 1.000000 

# Missing values propagate unless dropped.
c(kurtosis(c(y, NA)), kurtosis(c(y, NA), na.rm = TRUE))
#> [1]       NA 1.130882
```
