# Mean of the Skew Normal Distribution

Closed form: \\E\[Y\] = \mu + \sigma b \delta\\, with \\\delta =
\alpha/\sqrt{1+\alpha^2}\\ and \\b = \sqrt{2/\pi}\\. The location is not
the mean unless the shape is zero: a positive shape pulls the mass to
the right and the mean with it, by at most \\\sigma\sqrt{2/\pi} \approx
0.7979\sigma\\.

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

A numeric vector of means, of length equal to the longest of the three
components.

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale and \\\alpha\\ the
shape, in Azzalini's direct parametrization.

## See also

[`variance.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.SkewNormal1Distrib.md),
[`skewness.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.SkewNormal1Distrib.md),
[`skewnormal_delta()`](https://statmodels7.github.io/distributions7/reference/skewnormal_delta.md)
for the shared factor,
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md).

## Examples

``` r
d <- skewnormal1_distrib()

# The location plus sigma b delta.
delta <- 3 / sqrt(1 + 9)
all.equal(mean(d, list(mu = 0, sigma = 1, alpha = 3)),
          sqrt(2 / pi) * delta)
#> [1] TRUE

# At shape zero the family is Gaussian and the location is the mean.
mean(d, list(mu = 2, sigma = 1, alpha = 0))
#> [1] 2
```
