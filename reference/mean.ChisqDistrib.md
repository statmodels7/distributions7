# Mean of the Chi-Squared Distribution

Closed form: \\E\[Y\] = \mu\\. The family carries one parameter, which
is both the mean and the degrees of freedom, the two coinciding for this
law. It is not restricted to whole numbers here, so the family is the
gamma with variance tied to twice the mean.

## Arguments

- x:

  A `ChisqDistrib`, from
  [`chisq_distrib()`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md).

- theta:

  A named list with one component, `mu` (the mean and the degrees of
  freedom, positive), a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, the length of `theta$mu`.

## See also

[`variance.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.ChisqDistrib.md),
which is twice this;
[`mean.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Gamma2Distrib.md)
for the containing family;
[`chisq_distrib()`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md).

## Examples

``` r
d <- chisq_distrib()

# The one parameter is the mean and the degrees of freedom at once.
mean(d, list(mu = c(1, 5, 20)))
#> [1]  1  5 20
```
