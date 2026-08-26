# Multinomial Marginal

Returns the law of one coordinate, which is \\\mathrm{Binomial}(n,
p_j)\\: the other categories collapse into a single failure, and the
trial count is unchanged. The returned object carries the same `size` as
the multinomial and its `mu` is \\p_j\\.

Several coordinates at once are refused. A sub-vector is again
multinomial, but only after the remaining outcomes are collapsed into a
category of their own, and returning that object under this name would
mislead. The error says so.

## Arguments

- distrib:

  A `MultinomialDistrib` object, from
  [`multinomial_distrib()`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md).

- theta:

  A named list of the simplex's free values on the parameter scale, each
  of length 1.

- which:

  A single integer in \\1, \dots, p\\, the category wanted. A vector of
  length other than 1 signals an error explaining why.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with `distrib`, a
[`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md)
object at the same `size`, and `theta`, a list holding the marginal's
`mu`.

## See also

[`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md)
for the family returned,
[`mv_sigma.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.MultinomialDistrib.md)
for the coordinate variances, and
[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md)
for the generic.

## Examples

``` r
d <- multinomial_distrib(3, size = 5)
th <- list(probs_alr1 = 0.3, probs_alr2 = -0.2)
m <- mv_marginal(d, th, which = 2)
m$distrib
#> Distribution: Binomial
#> Type:         Discrete
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu (probability)        | Link: logit      | Domain: (0, 1)
m$theta
#> $mu
#> [1] 0.2583897
#> 

# The marginal's mean and variance are the second coordinate's.
rbind(marginal = c(mean(m$distrib, m$theta), variance(m$distrib, m$theta)),
      joint = c(mv_location(d, th)[2], mv_sigma(d, th)[2, 2]))
#>              [,1]      [,2]
#> marginal 1.291948 0.9581222
#> joint    1.291948 0.9581222

# Every marginal carries the same trial count.
vapply(1:3, function(j) mv_marginal(d, th, which = j)$distrib@size,
       numeric(1))
#> [1] 5 5 5
```
