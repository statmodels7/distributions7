# Skewness of the Beta Distribution in Two Shapes

Closed form: \$\$\gamma_1 = \frac{2(\beta - \alpha)\sqrt{\alpha +
\beta + 1}} {(\alpha + \beta + 2)\sqrt{\alpha\beta}}.\$\$ It takes the
sign of \\\beta - \alpha\\, so a beta is right-skewed when the first
shape is the smaller and exactly symmetric when the two are equal.

## Arguments

- x:

  A `Beta2Distrib`, from
  [`beta2_distrib()`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md).

- theta:

  A named list with components `alpha` and `beta`, both positive, each a
  numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector, of length
`max(length(theta$alpha), length(theta$beta))`.

## Notation

\\\alpha \> 0\\ and \\\beta \> 0\\ are the two shapes.

## See also

[`kurtosis.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Beta2Distrib.md),
written in the same two shapes;
[`skewness.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Beta1Distrib.md);
[`beta2_distrib()`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md).

## Examples

``` r
d <- beta2_distrib()

# Zero when the two shapes are equal.
skewness(d, list(alpha = 3, beta = 3))
#> [1] 0

# It agrees with the mean-precision parametrization on the same law.
all.equal(skewness(d, list(alpha = 2, beta = 3)),
          skewness(beta1_distrib(), list(mu = 0.4, phi = 5)))
#> [1] TRUE
```
