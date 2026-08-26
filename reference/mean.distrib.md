# Mean of a Distribution Object

Evaluates \\E\[Y\]\\ by quadrature over the support of a continuous
family or by summing the series of a discrete one, through one call of
[`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md)
at `p = 1`. This is the fallback that runs when the family has
registered no closed form. Of the 45 shipped families only the two von
Mises reach it; the rest carry a formula, which is a hundredfold cheaper
and exact.

## Arguments

- x:

  An object inheriting from `distrib`.

- theta:

  A named list of parameters on the parameter scale, one component per
  parameter of `x`. Components may be vectors, in which case one mean is
  returned per setting.

- ...:

  Passed to
  [`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md)
  and from there to
  [`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md).

## Value

A numeric vector of means, of length equal to the longest component of
`theta`. `NaN` or `Inf` where the integral does not converge, with no
warning: a family whose mean does not exist registers its own method
instead of relying on this one.

## See also

[`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md)
for the quadrature,
[`mean.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.CauchyDistrib.md)
for a family that overrides this method because its mean does not exist,
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md),
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md)
and
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md)
for the higher moments.

## Examples

``` r
# The von Mises is one of the two families with no closed-form mean, so this
# method is what answers for it. On (-pi, pi] with mu = 0 the mean is 0.
round(mean(vonmises1_distrib(), list(mu = 0, kappa = 2)), 12)
#> [1] 0

# A family with a closed form never reaches here; the numbers agree anyway.
all.equal(moment(gaussian1_distrib(), list(mu = 2, sigma = 3), p = 1), 2)
#> [1] TRUE
```
