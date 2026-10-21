# Sample Variance

Computes the sample variance of a numeric vector by delegating to
[`stats::var()`](https://rdrr.io/r/stats/cor.html), which divides by
\\n - 1\\. The method exists so that
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
reads the same on data as on a distribution object; the method computes
only what [`stats::var()`](https://rdrr.io/r/stats/cor.html) computes.

## Arguments

- x:

  A numeric vector. A vector of length 1 gives `NA`, as
  [`stats::var()`](https://rdrr.io/r/stats/cor.html) does, because zero
  degrees of freedom remain.

- na.rm:

  Should missing values be dropped before the variance is taken? It is a
  single logical, `FALSE` by default, and `FALSE` propagates `NA`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A single number, `NA` if `x` has fewer than two non-missing values.

## See also

[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
for the generic,
[`std_dev.numeric()`](https://statmodels7.github.io/distributions7/reference/std_dev.numeric.md)
for the square root,
[`skewness.numeric()`](https://statmodels7.github.io/distributions7/reference/skewness.numeric.md)
and
[`kurtosis.numeric()`](https://statmodels7.github.io/distributions7/reference/kurtosis.numeric.md),
which use the \\n\\ denominator instead.

## Examples

``` r
set.seed(1)
y <- rnorm(50)
all.equal(variance(y), var(y))
#> [1] TRUE

# The denominator is n - 1, not n.
all.equal(variance(y), sum((y - mean(y))^2) / (length(y) - 1))
#> [1] TRUE

# Missing values propagate unless dropped.
c(variance(c(y, NA)), variance(c(y, NA), na.rm = TRUE))
#> [1]        NA 0.6912159
```
