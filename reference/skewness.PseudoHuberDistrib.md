# Skewness of the Pseudo-Huber Distribution

Exactly zero at every parameter value. The density is symmetric about
\\\mu\\, so every odd central moment vanishes and the third standardized
one with it. The constant is returned directly, recycled to the length
the parameters imply, so no quadrature is run.

## Arguments

- x:

  A `PseudoHuberDistrib`, from
  [`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md).

- theta:

  A named list with components `mu`, `sigma` and `nu`, each a numeric
  vector of length 1 or `n`. The values are not read, only their
  lengths, but the list is still aligned and validated.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of zeros, of length equal to the longest of the three
components.

## See also

[`kurtosis.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.PseudoHuberDistrib.md),
which is not constant;
[`mean.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.PseudoHuberDistrib.md);
[`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md).

## Examples

``` r
d <- pseudohuber_distrib()

# Zero at every shape, by symmetry, with one value per setting.
skewness(d, list(mu = 0, sigma = 1, nu = c(1, 2, 3)))
#> [1] 0 0 0
```
