# Binomial Random Number Generator

Draws `n` independent counts of successes by calling
[`stats::rbinom()`](https://rdrr.io/r/stats/Binomial.html) at
`size = distrib@size`, so the draws come from R's own generator and
depend on `.Random.seed` in the usual way.

## Arguments

- distrib:

  A `BinomialDistrib` object, from
  [`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md).
  Its `size` property supplies the number of trials, and may be one
  value per draw.

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of length `n`. A value of length 1 is recycled. `mu` must
  lie in \\(0, 1)\\.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

An integer vector of `n` counts, each between 0 and the corresponding
`size`.

## See also

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the probability back, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- binomial_distrib(size = 10)

# Same generator as stats::rbinom, so the same seed gives the same draws.
set.seed(2)
a <- distrib_rng(d, 5, list(mu = 0.3))
set.seed(2)
identical(a, rbinom(5, size = 10, prob = 0.3))
#> [1] TRUE

# The sample mean estimates n p and the variance n p (1-p).
set.seed(4)
z <- distrib_rng(d, 2e4, list(mu = 0.3))
c(mean = mean(z), n_p = 10 * 0.3,
  var = var(z), n_p_q = 10 * 0.3 * 0.7)
#>     mean      n_p      var    n_p_q 
#> 2.991200 3.000000 2.089827 2.100000 

# Unequal group sizes, one per draw.
set.seed(5)
g <- binomial_distrib(size = c(5, 10, 20, 50))
rbind(size = c(5, 10, 20, 50), draw = distrib_rng(g, 4, list(mu = 0.3)))
#>      [,1] [,2] [,3] [,4]
#> size    5   10   20   50
#> draw    1    4    9   13
```
