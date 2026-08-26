# CDF Gradient When Only Some Parameters Are Location-Scale

The
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
body for a family that is location-scale in \\(\mu, \sigma)\\ and
carries a further shape parameter: the two location-scale directions in
closed form, \\-f\\ and \\-z f\\, and the shape directions from
[`numerical_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv.md).
The Student t and the pseudo-Huber are the two families that register
it.

## Usage

``` r
partial_loc_scale_grad_cdf(
  distrib,
  q,
  theta,
  lower.tail = TRUE,
  log = TRUE,
  ...
)
```

## Arguments

- distrib:

  An object inheriting from `distrib` whose first two parameters are a
  location and a scale, and which has at least one more.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters, the location first and the scale second.

- lower.tail:

  Is the lower tail wanted? A single logical, `TRUE` by default.

- log:

  Are derivatives of the log probability wanted? A single logical,
  `TRUE` by default.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors, one per parameter, in the parameter
order.

## Why the shape is differenced

In both families the shape direction is the derivative of a
hypergeometric-type integral with respect to its parameter, and has no
elementary form. Only those components are passed to the numerical
route, through its `which` argument, so the differencing costs cdf
evaluations for the shape alone.

## What the split is worth

It is a speed gain, and the size of it depends on how dear the family's
cdf is. On a pseudo-Huber, whose distribution function is itself a
quadrature, a gradient at 500 quantiles costs 0.08 s here against 0.18 s
when every component is differenced. On accuracy the two agree closely:
measured over shapes from 0.1 to 100 and quantiles out to eight scales,
the differenced mean component is within \\2\times10^{-11}\\ to
\\3\times10^{-8}\\ relative of the closed one.

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale, \\z =
(q-\mu)/\sigma\\ and \\f\\ the density.

## See also

[`loc_scale_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_grad_cdf.md),
the body for a family with no shape;
[`numerical_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv.md)
for the shape components;
[`distrib_grad_cdf.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.StudentT1Distrib.md)
and
[`distrib_grad_cdf.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.PseudoHuberDistrib.md),
the two registrations.
