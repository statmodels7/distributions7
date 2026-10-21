# Excess Kurtosis of the Beta Distribution in Two Shapes

The excess kurtosis has the closed form \$\$\gamma_2 =
\frac{6\\(\alpha-\beta)^2(\alpha+\beta+1) -
\alpha\beta(\alpha+\beta+2)\\}
{\alpha\beta(\alpha+\beta+2)(\alpha+\beta+3)}.\$\$ It is negative over a
wide range of shapes, because the support is bounded at both ends: at
\\\alpha = \beta = 1\\ the density is uniform and the value is \\-6/5\\,
and as both shapes go to zero it tends to \\-2\\, the smallest excess
kurtosis that any distribution can have. The method evaluates the form
in the mean and the precision of
[`beta_shape_moments()`](https://statmodels7.github.io/distributions7/reference/beta_shape_moments.md),
in which \\\alpha\beta\\ is not formed.

## Arguments

- x:

  A `Beta2Distrib`, from
  [`beta2_distrib()`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md).

- theta:

  A named list with components `alpha` and `beta`, both positive; each
  is a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of excess kurtoses, of length
`max(length(theta$alpha), length(theta$beta))`.

## Notation

\\\alpha \> 0\\ and \\\beta \> 0\\ are the two shapes.

## See also

[`skewness.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Beta2Distrib.md),
written in the same two shapes;
[`kurtosis.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Beta1Distrib.md);
[`beta2_distrib()`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md).

## Examples

``` r
d <- beta2_distrib()

# The uniform case has excess kurtosis -6/5.
all.equal(kurtosis(d, list(alpha = 1, beta = 1)), -1.2)
#> [1] TRUE

# It agrees with the mean-precision parametrization on the same law.
all.equal(kurtosis(d, list(alpha = 2, beta = 3)),
          kurtosis(beta1_distrib(), list(mu = 0.4, phi = 5)))
#> [1] TRUE
```
