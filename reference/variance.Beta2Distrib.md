# Variance of the Beta Distribution in Two Shapes

Closed form: \$\$\operatorname{Var}(Y) = \frac{\alpha\beta}
{(\alpha+\beta)^2(\alpha+\beta+1)}.\$\$ Writing \\\mu =
\alpha/(\alpha+\beta)\\ and \\\phi = \alpha+\beta\\ this is the
\\\mu(1-\mu)/(\phi+1)\\ of
[`variance.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Beta1Distrib.md),
so the shapes and the mean-precision pair describe the same law.

## Arguments

- x:

  A `Beta2Distrib`, from
  [`beta2_distrib()`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md).

- theta:

  A named list with components `alpha` and `beta`, both positive, each a
  numeric vector of length 1 or `n`. Both shapes below 1 puts the mass
  at the two endpoints.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of length
`max(length(theta$alpha), length(theta$beta))`.

## Notation

\\\alpha \> 0\\ and \\\beta \> 0\\ are the two shapes.

## See also

[`mean.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Beta2Distrib.md),
[`skewness.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Beta2Distrib.md);
[`variance.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Beta1Distrib.md);
[`beta2_distrib()`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md).

## Examples

``` r
d <- beta2_distrib()

# The published form, written out.
all.equal(variance(d, list(alpha = 2, beta = 3)), 2 * 3 / (25 * 6))
#> [1] TRUE

# Both shapes 1 is the uniform, of variance 1/12.
all.equal(variance(d, list(alpha = 1, beta = 1)), 1 / 12)
#> [1] TRUE
```
