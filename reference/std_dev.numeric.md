# Sample Standard Deviation

The sample standard deviation of a numeric vector, delegated to
[`stats::sd()`](https://rdrr.io/r/stats/sd.html) and so the root of the
\\n - 1\\ variance. The method exists so that
[`std_dev()`](https://statmodels7.github.io/distributions7/reference/std_dev.md)
reads the same on data as on a distribution object.

## Arguments

- x:

  A numeric vector. A vector of length 1 gives `NA`, as
  [`stats::sd()`](https://rdrr.io/r/stats/sd.html) does.

- na.rm:

  Should missing values be dropped before the standard deviation is
  taken? A single logical, `FALSE` by default, which propagates `NA`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A single number, `NA` if `x` has fewer than two non-missing values.

## See also

[`std_dev()`](https://statmodels7.github.io/distributions7/reference/std_dev.md)
for the generic,
[`variance.numeric()`](https://statmodels7.github.io/distributions7/reference/variance.numeric.md)
for the square,
[`skewness.numeric()`](https://statmodels7.github.io/distributions7/reference/skewness.numeric.md)
and
[`kurtosis.numeric()`](https://statmodels7.github.io/distributions7/reference/kurtosis.numeric.md).

## Examples

``` r
set.seed(1)
y <- rnorm(50)
all.equal(std_dev(y), sd(y))
#> [1] TRUE
all.equal(std_dev(y), sqrt(variance(y)))
#> [1] TRUE

# Missing values propagate unless dropped.
c(std_dev(c(y, NA)), std_dev(c(y, NA), na.rm = TRUE))
#> [1]        NA 0.8313939
```
