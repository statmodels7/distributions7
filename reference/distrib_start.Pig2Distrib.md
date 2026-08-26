# Orthogonal Poisson-Inverse Gaussian Starting Values

Returns the moment estimate of
[`pig1()`](https://statmodels7.github.io/distributions7/reference/distrib_start.Pig1Distrib.md)
carried onto this chart. The two families are the same law, so the
estimate is the same estimate: \\\hat\mu = \bar y\\ and \\\hat\sigma =
(s^2 - \bar y)/\bar y^2\\ as before, then \$\$\alpha = \frac{\sqrt{1 +
2\sigma\mu}}{\sigma}.\$\$

The same two floors apply, so an underdispersed sample gives a very
large \\\alpha\\, which is where this parametrization puts the Poisson
limit.

## Arguments

- distrib:

  A `Pig2Distrib` object.

- y:

  A numeric vector of counts.

- n_start:

  Ignored: one moment start is returned.

- ...:

  Unused.

## Value

A list of length 1 holding one named parameter list with components `mu`
and `alpha`.

## See also

[`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md)
for the family;
[`distrib_start.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_start.Pig1Distrib.md)
for the estimate this maps;
[`distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.md)
for the generic.

## Examples

``` r
set.seed(4)
d <- pig2_distrib()
alpha0 <- sqrt(1 + 2 * 0.5 * 4) / 0.5
y <- distrib_rng(d, 5000, list(mu = 4, alpha = alpha0))

# The inversion recovers what the sample was drawn from.
rbind(start = unlist(distrib_start(d, y)[[1]]),
      truth = c(mu = 4, alpha = alpha0))
#>           mu    alpha
#> start 3.9834 4.528028
#> truth 4.0000 4.472136

# It is pig1's estimate carried through the map.
s <- (var(y) - mean(y)) / mean(y)^2
c(alpha = sqrt(1 + 2 * s * mean(y)) / s)
#>    alpha 
#> 4.528028 
```
