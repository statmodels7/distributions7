# Generalized Gamma Random Generation

Draws by the family's own defining representation: a Gamma variate with
shape \\d/p\\ and unit rate, raised to the power \\1/p\\ and multiplied
by \\a\\. It is exact and costs one `rgamma` call, so neither inversion
nor the base class's ratio-of-uniforms fallback is involved.

## Arguments

- distrib:

  A `GenGamma1Distrib` object, from
  [`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with components `a`, `d` and `p`, each a numeric vector
  of length 1 or of length `n`; a component of length 1 is recycled, so
  a parameter varying by observation draws one value per observation
  from its own member of the family.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` strictly positive draws.

## See also

[`distrib_pdf.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.GenGamma1Distrib.md)
for the density the draws follow,
[`mean.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.GenGamma1Distrib.md)
for the moments they reproduce, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- gengamma1_distrib()
th <- list(a = 2, d = 3, p = 1.5)

set.seed(51)
x <- distrib_rng(d, 2e5, th)

# Two moments against their closed forms.
rbind(sample = c(mean(x), var(x)),
      theory = c(mean(d, th), variance(d, th)))
#>            [,1]     [,2]
#> sample 3.007083 2.054553
#> theory 3.009151 2.057644

# The representation, written out at the same seed.
set.seed(51)
all.equal(x, 2 * rgamma(2e5, shape = 3 / 1.5)^(1 / 1.5))
#> [1] TRUE
```
