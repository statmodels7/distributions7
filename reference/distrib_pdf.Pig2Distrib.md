# Orthogonal Poisson-Inverse Gaussian Probability Mass Function

Computes the same mass as
[`distrib_pdf.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Pig1Distrib.md),
at the dispersion \\\sigma = (\mu + \sqrt{\mu^2 + \alpha^2})/\alpha^2\\
that \\\alpha\\ implies. The parameter \\\alpha\\ of this
parametrization **is** the argument \\\sqrt{1 + 2\sigma\mu}/\sigma\\ at
which the Bessel function is evaluated, so in this coordinate the mass
function's own argument is a parameter.

The compiled kernel takes \\\alpha\\ directly, so
[`pig2_sigma()`](https://statmodels7.github.io/distributions7/reference/pig2_sigma.md)
is not called here and no map is composed.

## Arguments

- distrib:

  A `Pig2Distrib` object, from
  [`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md).

- y:

  A numeric vector of counts. A negative, non-integer or non-finite
  value is off the support and gives a probability of 0, or `-Inf` with
  `log = TRUE`.

- theta:

  A named list with components `mu` and `alpha`, each a numeric vector
  of length 1 or of the length of `y`. Both must be strictly positive.

- log:

  Logical of length 1. When `TRUE` the log-probability is returned.
  Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the compiled kernel may
  use. Defaults to `1L`. The result does not depend on the count.

## Value

A numeric vector of probabilities, of the length of the recycled inputs.

## Notation

\\\mu\\ is the mean, \\\alpha\\ the Bessel argument, and \\\sigma\\ the
dispersion of
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md),
which \\\alpha\\ determines.

## See also

[`distrib_pdf.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Pig1Distrib.md)
for the same mass in mean and dispersion,
[`pig2_sigma()`](https://statmodels7.github.io/distributions7/reference/pig2_sigma.md)
for the map, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- pig2_distrib()
d1 <- pig1_distrib()
y <- 0:6
al <- sqrt(1 / 0.8^2 + 2 * 3 / 0.8)

# The same law, reached through the map.
all.equal(distrib_pdf(d, y, list(mu = 3, alpha = al)),
          distrib_pdf(d1, y, list(mu = 3, sigma = 0.8)))
#> [1] TRUE

# The mass sums to one.
sum(distrib_pdf(d, 0:300, list(mu = 3, alpha = al)))
#> [1] 1

# A large alpha is a small dispersion: the mass tends to the Poisson's.
rbind(pig2 = distrib_pdf(d, y, list(mu = 3, alpha = 1e4)),
      poisson = dpois(y, 3))
#>               [,1]      [,2]      [,3]      [,4]      [,5]      [,6]       [,7]
#> pig2    0.04980948 0.1493836 0.2240306 0.2240082 0.1680061 0.1008138 0.05041696
#> poisson 0.04978707 0.1493612 0.2240418 0.2240418 0.1680314 0.1008188 0.05040941
```
