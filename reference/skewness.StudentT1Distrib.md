# Skewness of the Student t Distribution

The skewness is zero for \\\nu \> 3\\ and `NaN` at or below three, where
the third moment does not exist. The density is symmetric about \\\mu\\,
so the value is zero wherever it is defined.

## Arguments

- x:

  A `StudentT1Distrib`, from
  [`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

- theta:

  A named list with components `mu`, `sigma` (positive) and `nu`
  (positive), in which each component is a numeric vector of length 1 or
  `n`. Only `nu` enters the value; settings with \\\nu \le 3\\ give
  `NaN`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector, whose length is the one that the recycled parameters
imply, with zero where \\\nu \> 3\\ and `NaN` elsewhere.

## Notation

\\\mu\\ is the location and \\\nu \> 0\\ the degrees of freedom.

## See also

[`kurtosis.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.StudentT1Distrib.md),
whose threshold is 4;
[`skewness.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.SkewTDistrib.md),
the asymmetric extension;
[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

## Examples

``` r
d <- student_t1_distrib()

# The skewness is undefined at or below three degrees of freedom, and zero by
# symmetry above.
skewness(d, list(mu = 0, sigma = 1, nu = c(2, 3, 4)))
#> [1] NaN NaN   0
```
