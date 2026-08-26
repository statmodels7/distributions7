# Transformed Probability Density Function

Computes the change-of-variables density \$\$f_Y(y) =
f_X(g^{-1}(y))\\\lvert J(y)\rvert,\$\$ entirely on the LOG scale and
exponentiated at the end, so that a density and a Jacobian which would
each overflow or underflow on their own combine as a sum of
representable numbers.

## Arguments

- distrib:

  A `TransformedDistrib` object, from
  [`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md).

- y:

  A numeric vector of observations on the TRANSFORMED scale. A point
  outside the transformed support gives `0`, the parent's density at its
  preimage being zero there.

- theta:

  A named list of the PARENT's parameters, unchanged by the
  transformation.

- log:

  Logical of length 1. When `TRUE` the log-density is returned, which is
  the quantity actually computed. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of the recycled length of `y` and `theta`.

## Details

Two singular cases are handled explicitly. A parent log-density that is
`+Inf`, which a density with a pole produces, is clamped to
`log(.Machine$double.xmax)`, so that a quadrature over the transformed
density meets a large number instead of `Inf - Inf`. A parent density of
exactly zero WINS over an infinite Jacobian, so the result is `-Inf`.
The point is outside the transformed support and `-Inf` says so, where
the unguarded sum would give `NaN`.

## Notation

\\g\\ is the transformation, \\J\\ the Jacobian of its inverse and
\\f_X\\, \\f_Y\\ the parent's and the transformed density.

## See also

[`distrib_cdf.TransformedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.TransformedDistrib.md)
for the distribution function,
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md)
for the family, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- transformation(gaussian1_distrib(), exp_transform())
theta <- list(mu = 0.5, sigma = 0.8)
y <- c(0.5, 1, 3)

distrib_pdf(d, y, theta)
#> [1] 0.3279692 0.4102012 0.1256371

# The exponential of a gaussian is the lognormal.
all.equal(distrib_pdf(d, y, theta), dlnorm(y, 0.5, 0.8))
#> [1] TRUE

# It is a density: it integrates to one over the transformed support.
integrate(function(z) distrib_pdf(d, z, theta), 0, Inf)$value
#> [1] 1

# Written out, the parent's density at the preimage times the Jacobian.
tr <- d@transformer
all.equal(distrib_pdf(d, y, theta),
          dnorm(tr@trans_inv(y), 0.5, 0.8) * tr@trans_abs_jac(y))
#> [1] "Mean relative difference: 1.347779"
```
