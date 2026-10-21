# Mean of the Beta-Binomial Distribution in Two Shapes

The mean has the closed form \\E\[Y\] = n\alpha/(\alpha+\beta)\\, the
trial count times the mean of the mixing beta. This is the classical
two-shape parametrization;
[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md)
carries the same law with \\\mu = \alpha/(\alpha+\beta)\\ and \\\sigma =
1/(\alpha+\beta)\\.

## Arguments

- x:

  A `BetaBinom2Distrib`, from
  [`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md),
  carrying the trial count in its `size` property.

- theta:

  A named list with components `alpha` and `beta`, both positive; each
  is a numeric vector of length 1 or of the common length of the
  parameter settings.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, of length
`max(length(theta$alpha), length(theta$beta))`.

## Notation

\\\alpha \> 0\\ and \\\beta \> 0\\ are the shapes of the mixing beta and
\\n\\ the number of trials, held on the object.

## See also

[`variance.BetaBinom2Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.BetaBinom2Distrib.md),
where the dispersion does enter;
[`mean.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.BetaBinom1Distrib.md)
for the mean-dispersion parametrization;
[`betabinom_moments()`](https://statmodels7.github.io/distributions7/reference/betabinom_moments.md);
[`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md).

## Examples

``` r
d <- betabinom2_distrib(size = 10)

# The mean is ten trials times the mixing beta's mean of 0.4.
all.equal(mean(d, list(alpha = 2, beta = 3)), 4)
#> [1] TRUE
```
