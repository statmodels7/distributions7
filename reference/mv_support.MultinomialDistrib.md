# The Support Points of a Multinomial

Enumerates the whole support: every vector of \\p\\ non-negative
integers summing to the trial count \\n\\, one row per point. These are
the weak compositions of \\n\\ into \\p\\ parts, supplied by
[`numericals7::compositions()`](https://statmodels7.github.io/numericals7/reference/compositions.html).

Because the support is enumerated, every expectation is an exact sum.
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
uses it to test that the mass adds to one, an equality with tolerance
`1e-10`; an importance-sampling check compares only against Monte Carlo
error, and a normalization wrong by a thousandth would pass it.

## Arguments

- distrib:

  A `MultinomialDistrib` object, from
  [`multinomial_distrib()`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md),
  read for its `size` and `n_dim`.

- theta:

  Ignored. The support depends only on \\n\\ and \\p\\, so any value,
  including `NULL`, gives the same answer.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric matrix with \\\binom{n+p-1}{p-1}\\ rows and \\p\\ columns. The
count grows quickly in both arguments: 21 rows at \\n = 5, p = 3\\, and
1001 at \\n = 10, p = 5\\.

## See also

[`numericals7::compositions()`](https://statmodels7.github.io/numericals7/reference/compositions.html)
for the enumeration,
[`distrib_pdf.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.MultinomialDistrib.md)
for the mass on these points,
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md),
which consumes this, and
[`mv_support()`](https://statmodels7.github.io/distributions7/reference/mv_support.md)
for the generic.

## Examples

``` r
d <- multinomial_distrib(3, size = 5)
supp <- mv_support(d, NULL)
dim(supp)
#> [1] 21  3
head(supp)
#>      [,1] [,2] [,3]
#> [1,]    0    0    5
#> [2,]    0    1    4
#> [3,]    0    2    3
#> [4,]    0    3    2
#> [5,]    0    4    1
#> [6,]    0    5    0

# Every row sums to the trial count, and the count of rows is the number of
# weak compositions.
c(all_sum_to_5 = all(rowSums(supp) == 5),
  rows = nrow(supp), formula = choose(5 + 3 - 1, 3 - 1))
#> all_sum_to_5         rows      formula 
#>            1           21           21 

# The mass over the support adds to one exactly.
sum(distrib_pdf(d, supp, list(probs_alr1 = 0.3, probs_alr2 = -0.2))) - 1
#> [1] -7.771561e-16
```
