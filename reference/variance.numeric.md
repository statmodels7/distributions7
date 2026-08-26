# Sample Variance

The sample variance of a numeric vector, delegated to
[`stats::var()`](https://rdrr.io/r/stats/cor.html) and so divided by
\\n - 1\\. The method exists so that
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
reads the same on data as on a distribution object; nothing is computed
here that [`stats::var()`](https://rdrr.io/r/stats/cor.html) does not
compute.

## Arguments

- x:

  A numeric vector. A vector of length 1 gives `NA`, as
  [`stats::var()`](https://rdrr.io/r/stats/cor.html) does, there being
  no degree of freedom left.

- na.rm:

  Should missing values be dropped before the variance is taken? A
  single logical, `FALSE` by default, which propagates `NA`.

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

# The n - 1 denominator, not n.
all.equal(variance(y), sum((y - mean(y))^2) / (length(y) - 1))
#> [1] TRUE

# Missing values propagate unless dropped.
c(variance(c(y, NA)), variance(c(y, NA), na.rm = TRUE))
#> [1]        NA 0.6912159
```
