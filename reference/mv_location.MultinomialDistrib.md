# Multinomial Mean Vector

Returns \\\mathbb{E}\[Y\] = n p\\, the trial count times the probability
vector, so the result sums to \\n\\ and not to one. Divide by the
object's `size` to recover the probabilities themselves.
[`mean.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.MultinomialDistrib.md)
delegates here.

## Arguments

- distrib:

  A `MultinomialDistrib` object, from
  [`multinomial_distrib()`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md).

- theta:

  A named list of the simplex's free values on the parameter scale, each
  of length 1.

## Value

A numeric vector of length \\p\\, strictly positive and summing to the
object's `size`.

## See also

[`mv_sigma.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.MultinomialDistrib.md)
for the covariance,
[`mean.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.MultinomialDistrib.md),
which calls this, and
[`mv_location()`](https://statmodels7.github.io/distributions7/reference/mv_location.md)
for the generic.

## Examples

``` r
d <- multinomial_distrib(3, size = 5)
th <- list(probs_alr1 = 0.3, probs_alr2 = -0.2)
m <- mv_location(d, th)
m
#> [1] 2.130063 1.291948 1.577989
sum(m)
#> [1] 5

# Dividing by the trial count gives the probabilities, which sum to one.
m / d@size
#> [1] 0.4260125 0.2583897 0.3155978
```
