# Poisson-Inverse Gaussian Random Generation

Draws from the family's own mixture representation, exactly: \\\lambda\\
from the inverse Gaussian with mean \\\mu\\ and shape \\\mu/\sigma\\,
then \\Y \mid \lambda\\ from the Poisson with that rate. The inverse
Gaussian's variance is then \\\sigma\mu^2\\, the value the mixing
construction requires, and the marginal variance comes out as \\\mu +
\sigma\mu^2\\.

The inverse Gaussian draw comes from
[`statmod::rinvgauss`](https://rdrr.io/pkg/statmod/man/invgauss.html),
so the whole generator is two vectorized calls and no rejection step.

## Arguments

- distrib:

  A `Pig1Distrib` object, from
  [`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of length `n`; a component of length 1 is recycled, so
  a parameter varying by observation draws one value per observation
  from its own member of the family.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

An integer-valued numeric vector of `n` draws.

## See also

[`distrib_pdf.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Pig1Distrib.md)
for the mass the draws follow,
[`distrib_rng.Pig2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.Pig2Distrib.md)
for the same sampler in orthogonal coordinates, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- pig1_distrib()
th <- list(mu = 3, sigma = 0.8)

set.seed(61)
x <- distrib_rng(d, 2e5, th)

# Two moments against their closed forms.
rbind(sample = c(mean(x), var(x)),
      theory = c(mean(d, th), variance(d, th)))
#>            [,1]    [,2]
#> sample 2.996995 10.1801
#> theory 3.000000 10.2000

# The empirical mass against the exact one, over the head of the support.
rbind(sample = as.numeric(table(factor(x, levels = 0:6))) / 2e5,
      exact = distrib_pdf(d, 0:6, th))
#>             [,1]      [,2]      [,3]      [,4]       [,5]       [,6]       [,7]
#> sample 0.1731100 0.2149400 0.1762400 0.1278850 0.09066000 0.06239500 0.04293000
#> exact  0.1719763 0.2142278 0.1777529 0.1289567 0.08968701 0.06196187 0.04309807
```
