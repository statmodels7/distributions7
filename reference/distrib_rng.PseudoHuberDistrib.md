# Pseudo-Huber Random Number Generator

Draws `n` independent variates as a normal variance mixture, \\Y = \mu +
\sqrt{W} Z\\ with \\Z\\ standard normal and \\W = \sigma^2 \sqrt{\nu}\\
X\\, where \\X\\ follows the generalized inverse Gaussian law with
density proportional to \\\exp\\-\sqrt{\nu}\\(x + 1/x)/2\\\\. \\X\\ is
drawn by the ratio-of-uniforms method with the mode shifted to the
origin (Hormann and Leydold, 2014), whose bounding box is set once per
distinct \\\nu\\; below \\\nu = 10^{-16}\\ the mixing law is its
exponential limit. The draws depend on `.Random.seed` in the usual way.

## Arguments

- distrib:

  A `PseudoHuberDistrib` object, from
  [`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with components `mu`, `sigma` and `nu`, each a numeric
  vector of length 1 or of length `n`. A component of length 1 is
  recycled, so a vector of length `n` draws one variate per parameter
  setting. `sigma` and `nu` must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` draws.

## References

Hormann, W. and Leydold, J. (2014). Generating generalized inverse
Gaussian random variates. *Statistics and Computing*, **24**(4),
547-557.

## See also

[`distrib_quantile.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.PseudoHuberDistrib.md)
for the quantile function,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters back from a sample, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- pseudohuber_distrib()
th <- list(mu = 0.4, sigma = 1.2, nu = 2)

# The moments of a sample sit where the sampling error puts them.
set.seed(6)
z <- distrib_rng(d, 1e5, th)
rbind(sample = c(mean(z), var(z)),
      theoretical = c(mean(d, th), variance(d, th)))
#>                  [,1]     [,2]
#> sample      0.3954567 4.463978
#> theoretical 0.4000000 4.429997
```
