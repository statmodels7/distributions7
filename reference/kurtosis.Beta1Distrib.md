# Excess Kurtosis of the Beta Distribution

The excess kurtosis has the closed form, in the two shapes \\a =
\mu\phi\\ and \\b = (1-\mu)\phi\\, \$\$\gamma_2
=\frac{6\\(a-b)^2(a+b+1) - ab(a+b+2)\\} {ab(a+b+2)(a+b+3)}.\$\$ The
excess kurtosis is negative over a wide range of parameters because the
support is bounded. At \\\mu = 1/2\\, \\\phi = 2\\ the density is
uniform and the value is \\-6/5\\; at \\\mu = 1/2\\ it tends to \\-2\\,
the smallest value that any distribution can have, as \\\phi\\ goes to
zero. The method evaluates the form in \\(\mu, \phi)\\ of
[`beta_shape_moments()`](https://statmodels7.github.io/distributions7/reference/beta_shape_moments.md),
in which \\ab\\ is not formed.

## Arguments

- x:

  A `Beta1Distrib`, from
  [`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md).

- theta:

  A named list with components `mu` (strictly between 0 and 1) and `phi`
  (positive); each is a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of excess kurtoses, of length
`max(length(theta$mu), length(theta$phi))`.

## Notation

\\\mu \in (0,1)\\ is the mean, \\\phi \> 0\\ the precision, and \\a =
\mu\phi\\, \\b = (1-\mu)\phi\\ the two shapes.

## See also

[`skewness.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Beta1Distrib.md),
written in the same two shapes;
[`variance.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Beta1Distrib.md);
[`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md).

## Examples

``` r
d <- beta1_distrib()

# The uniform case has excess kurtosis -6/5; a lower precision gives lower
# values, down to -2.
all.equal(kurtosis(d, list(mu = 0.5, phi = 2)), -1.2)
#> [1] TRUE

# The value is negative over the symmetric middle and positive as the mass
# piles at an end.
round(kurtosis(d, list(mu = c(0.5, 0.1, 0.02), phi = 5)), 4)
#> [1] -0.7500  3.8214 29.4774
```
