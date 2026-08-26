# Multinomial Random Generation

Draws `n` independent multinomial counts by calling
[`stats::rmultinom()`](https://rdrr.io/r/stats/Multinom.html) and
transposing, so that one **row** of the result is one observation. R
returns one column per draw; the package's convention throughout the
multivariate families is one row per observation, and the transpose is
the whole of the difference. The draws depend on `.Random.seed` in the
usual way.

## Arguments

- distrib:

  A `MultinomialDistrib` object, from
  [`multinomial_distrib()`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md).

- n:

  A single positive integer, the number of draws. Note that the number
  of **trials** is the object's `size` property, not this argument.

- theta:

  A named list of the simplex's free values on the parameter scale, each
  of length 1.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric matrix with `n` rows and \\p\\ columns, each row a vector of
non-negative integers summing to the object's `size`.

## See also

[`distrib_pdf.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.MultinomialDistrib.md)
for the mass,
[`mv_support.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_support.MultinomialDistrib.md)
for the points it can land on,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the probabilities back from a sample, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- multinomial_distrib(3, size = 5)
th <- list(probs_alr1 = 0.3, probs_alr2 = -0.2)

# One row per observation, each summing to the trial count.
set.seed(1)
Y <- distrib_rng(d, 4, th)
Y
#>      [,1] [,2] [,3]
#> [1,]    1    1    3
#> [2,]    2    2    1
#> [3,]    1    3    1
#> [4,]    4    1    0
rowSums(Y)
#> [1] 5 5 5 5

# The sample recovers the mean vector and the coordinate variances.
set.seed(4)
Z <- distrib_rng(d, 2e4, th)
rbind(sample = c(colMeans(Z), diag(var(Z))),
      theoretical = c(mv_location(d, th), diag(mv_sigma(d, th))))
#>                 [,1]     [,2]     [,3]     [,4]      [,5]     [,6]
#> sample      2.124750 1.284400 1.590850 1.234349 0.9476640 1.084701
#> theoretical 2.130063 1.291948 1.577989 1.222629 0.9581222 1.079979
```
