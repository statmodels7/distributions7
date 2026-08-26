# Orthogonal Poisson-Inverse Gaussian Random Generation

Draws from the mixture representation of
[`distrib_rng.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.Pig1Distrib.md),
at the dispersion
[`pig2_sigma()`](https://statmodels7.github.io/distributions7/reference/pig2_sigma.md)
implies: \\\lambda\\ from the inverse Gaussian with mean \\\mu\\ and
shape \\\mu/\sigma\\, then \\Y \mid \lambda\\ from the Poisson. This is
the one method of the family that composes the map; the derivative
methods take \\\alpha\\ directly.

## Arguments

- distrib:

  A `Pig2Distrib` object, from
  [`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with components `mu` and `alpha`, each a numeric vector
  of length 1 or of length `n`; a component of length 1 is recycled, so
  a parameter varying by observation draws one value per observation
  from its own member of the family.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

An integer-valued numeric vector of `n` draws.

## See also

[`distrib_rng.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.Pig1Distrib.md)
for the same sampler in mean and dispersion,
[`pig2_sigma()`](https://statmodels7.github.io/distributions7/reference/pig2_sigma.md)
for the map it composes, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- pig2_distrib()
al <- sqrt(1 / 0.8^2 + 2 * 3 / 0.8)
th <- list(mu = 3, alpha = al)

set.seed(64)
x <- distrib_rng(d, 2e5, th)

# Two moments against their closed forms.
rbind(sample = c(mean(x), var(x)),
      theory = c(mean(d, th), variance(d, th)))
#>           [,1]     [,2]
#> sample 2.99463 10.09935
#> theory 3.00000 10.20000

# The empirical mass against the exact one, over the head of the support.
rbind(sample = as.numeric(table(factor(x, levels = 0:6))) / 2e5,
      exact = distrib_pdf(d, 0:6, th))
#>             [,1]      [,2]      [,3]      [,4]       [,5]       [,6]       [,7]
#> sample 0.1720950 0.2148850 0.1774300 0.1277900 0.09031500 0.06149500 0.04420000
#> exact  0.1719763 0.2142278 0.1777529 0.1289567 0.08968701 0.06196187 0.04309807
```
