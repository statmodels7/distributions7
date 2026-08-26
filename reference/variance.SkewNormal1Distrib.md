# Variance of the Skew Normal Distribution

Closed form: \\\operatorname{Var}(Y) = \sigma^2(1 - b^2\delta^2)\\, with
\\\delta = \alpha/\sqrt{1+\alpha^2}\\ and \\b = \sqrt{2/\pi}\\. The
bracket is at most 1 and at least \\1 - 2/\pi \approx 0.3634\\, so
skewing the family narrows it: the scale is an upper bound on the
standard deviation and is attained only at shape zero.

## Arguments

- x:

  A `SkewNormal1Distrib`, from
  [`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md).

- theta:

  A named list with components `mu` (the location), `sigma` (the scale,
  positive) and `alpha` (the shape, any sign), each a numeric vector of
  length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of the length the recycled parameters
imply. The scale and `alpha` enter the value, so a setting that varies
`mu` alone repeats one number.

## Notation

\\\sigma \> 0\\ is the scale and \\\alpha\\ the shape, in Azzalini's
direct parametrization.

## See also

[`mean.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.SkewNormal1Distrib.md),
[`kurtosis.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.SkewNormal1Distrib.md),
[`skewnormal_delta()`](https://statmodels7.github.io/distributions7/reference/skewnormal_delta.md),
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md).

## Examples

``` r
d <- skewnormal1_distrib()

# The variance falls from sigma^2 towards sigma^2 (1 - 2/pi) as the shape grows.
round(variance(d, list(mu = 0, sigma = 1, alpha = c(0, 1, 3, 1e6))), 6)
#> [1] 1.000000 0.681690 0.427042 0.363380
1 - 2 / pi
#> [1] 0.3633802
```
