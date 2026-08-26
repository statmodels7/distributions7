# Skewness of the Beta Distribution

Closed form in the two shapes \\a = \mu\phi\\ and \\b = (1-\mu)\phi\\:
\$\$\gamma_1 = \frac{2(b - a)\sqrt{a + b + 1}} {(a + b +
2)\sqrt{ab}}.\$\$ It takes the sign of \\b - a\\, so a beta is
right-skewed below a mean of one half, left-skewed above it, and exactly
symmetric at \\\mu = 1/2\\ whatever the precision.

## Arguments

- x:

  A `Beta1Distrib`, from
  [`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md).

- theta:

  A named list with components `mu` (strictly between 0 and 1) and `phi`
  (positive), each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector, of length `max(length(theta$mu), length(theta$phi))`.

## Notation

\\\mu \in (0,1)\\ is the mean, \\\phi \> 0\\ the precision, and \\a =
\mu\phi\\, \\b = (1-\mu)\phi\\ the two shapes.

## See also

[`kurtosis.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Beta1Distrib.md),
written in the same two shapes;
[`variance.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Beta1Distrib.md);
[`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md).

## Examples

``` r
d <- beta1_distrib()

# The published form in the two shapes, written out.
a <- 0.3 * 5; b <- 0.7 * 5
all.equal(skewness(d, list(mu = 0.3, phi = 5)),
          2 * (b - a) * sqrt(a + b + 1) / ((a + b + 2) * sqrt(a * b)))
#> [1] TRUE

# Exactly zero at a mean of one half, at any precision.
skewness(d, list(mu = 0.5, phi = c(0.5, 2, 50)))
#> [1] 0 0 0
```
