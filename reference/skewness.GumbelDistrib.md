# Skewness of the Gumbel Distribution

Constant: \\\gamma_1 = 12\sqrt6\\\zeta(3)/\pi^3 \approx 1.1395\\, with
\\\zeta(3)\\ Apery's constant. The Gumbel is location-scale with no
shape parameter, so its standardized moments are numbers and not
functions: every Gumbel has this same right skew, whatever its location
and scale.

## Arguments

- x:

  A `GumbelDistrib`, from
  [`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md).

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or `n`. The values are not read, only their lengths.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector, every element \\12\sqrt6\\\zeta(3)/\pi^3\\, of length
`max(length(theta$mu), length(theta$sigma))`.

## Details

\\\zeta(3) = 1.2020569031595942854\\ is written out in the body. Base R
has no function that returns it, and one number does not justify a
dependency.

## Notation

\\\zeta\\ is the Riemann zeta function; \\\mu\\ and \\\sigma \> 0\\ are
the location and the scale, neither of which enters the value.

## See also

[`kurtosis.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.GumbelDistrib.md),
the other constant;
[`mean.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.GumbelDistrib.md)
and
[`variance.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.GumbelDistrib.md),
which do move with the parameters;
[`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md).

## Examples

``` r
d <- gumbel_distrib()

# One number for the whole family.
skewness(d, list(mu = c(0, 5, -3), sigma = c(1, 7, 0.2)))
#> [1] 1.139547 1.139547 1.139547

# Against the published constant.
all.equal(skewness(d, list(mu = 0, sigma = 1)),
          12 * sqrt(6) * 1.2020569031595942854 / pi^3)
#> [1] TRUE
```
