# Poisson-Inverse Gaussian Starting Values

Returns the method-of-moments estimate. The family has mean \\\mu\\ and
variance \\\mu + \sigma\mu^2\\, so setting both equal to the sample
gives \\\hat\mu = \bar y\\ and \\\hat\sigma = (s^2 - \bar y)/\bar y^2\\
directly, with no root to find.

Both are floored: \\\mu\\ just above zero, and \\\sigma\\ at \\10^{-3}\\
when the sample is **underdispersed** and the inversion returns a
negative number. A Poisson sample is the case that produces it, and
\\10^{-3}\\ is where this family is nearly Poisson, which is the right
place to start from there.

## Arguments

- distrib:

  A `Pig1Distrib` object.

- y:

  A numeric vector of counts.

- n_start:

  Ignored: one moment start is returned.

- ...:

  Unused.

## Value

A list of length 1 holding one named parameter list with components `mu`
and `sigma`.

## See also

[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
for the family;
[`distrib_start.Pig2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_start.Pig2Distrib.md),
the same estimate on the orthogonal chart;
[`distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.md)
for the generic.

## Examples

``` r
set.seed(4)
d <- pig1_distrib()
y <- distrib_rng(d, 5000, list(mu = 4, sigma = 0.5))

# The inversion recovers the parameters it was drawn from.
unlist(distrib_start(d, y)[[1]])
#>       mu    sigma 
#> 3.983400 0.488425 

# It is exactly the sample moments, read through mean and variance.
c(mu = mean(y), sigma = (var(y) - mean(y)) / mean(y)^2)
#>       mu    sigma 
#> 3.983400 0.488425 

# An underdispersed sample would give a negative sigma, so it is floored
# where the family is nearly Poisson.
set.seed(5)
unlist(distrib_start(d, rpois(2000, 4))[[1]])
#>    mu sigma 
#> 4.023 0.001 
```
