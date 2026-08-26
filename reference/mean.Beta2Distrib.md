# Mean of the Beta Distribution in Two Shapes

Closed form: \\E\[Y\] = \alpha/(\alpha+\beta)\\. This is the beta in its
classical two-shape parametrization, so the mean is a ratio rather than
a parameter;
[`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md)
carries it directly, with \\\mu = \alpha/(\alpha+\beta)\\ and \\\phi =
\alpha+\beta\\.

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

A numeric vector of means, of length
`max(length(theta$alpha), length(theta$beta))`, strictly inside
\\(0,1)\\.

## Notation

\\\alpha \> 0\\ and \\\beta \> 0\\ are the two shapes.

## See also

[`variance.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Beta2Distrib.md);
[`mean.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Beta1Distrib.md)
for the mean-precision parametrization;
[`beta2_distrib()`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md).

## Examples

``` r
d <- beta2_distrib()

# The ratio of the first shape to their sum.
all.equal(mean(d, list(alpha = 2, beta = 3)), 0.4)
#> [1] TRUE
```
