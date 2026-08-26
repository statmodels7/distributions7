# Bernoulli Random Number Generator

Draws `n` independent zero-one variates by calling
[`stats::rbinom()`](https://rdrr.io/r/stats/Binomial.html) at
`size = 1`, so the draws come from R's own generator and depend on
`.Random.seed` in the usual way.

## Arguments

- distrib:

  A `BernoulliDistrib` object, from
  [`bernoulli_distrib()`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of length `n`. A value of length 1 is recycled, so a
  vector of length `n` draws one variate per probability, the shape a
  regression on a Bernoulli response supplies. `mu` must lie in \\(0,
  1)\\.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

An integer vector of `n` zeros and ones.

## See also

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the probability back, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- bernoulli_distrib()

# Same generator as stats::rbinom, so the same seed gives the same draws.
set.seed(2)
a <- distrib_rng(d, 5, list(mu = 0.3))
set.seed(2)
identical(a, rbinom(5, size = 1, prob = 0.3))
#> [1] TRUE

# The sample proportion estimates mu.
set.seed(4)
z <- distrib_rng(d, 2e4, list(mu = 0.3))
c(proportion = mean(z), variance = var(z), p_times_q = 0.3 * 0.7)
#> proportion   variance  p_times_q 
#>  0.3010500  0.2104294  0.2100000 

# A probability per observation, which is what a regression supplies.
set.seed(5)
p <- plogis(seq(-2, 2, length.out = 6))
rbind(p = round(p, 3), draw = distrib_rng(d, 6, list(mu = p)))
#>       [,1]  [,2]  [,3]  [,4]  [,5]  [,6]
#> p    0.119 0.231 0.401 0.599 0.769 0.881
#> draw 0.000 0.000 1.000 1.000 1.000 1.000
```
