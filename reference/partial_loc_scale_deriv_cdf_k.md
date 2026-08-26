# Higher Log-CDF Derivatives When Only Some Parameters Are Location-Scale

Builds the
[`distrib_deriv3_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_cdf.md)
or
[`distrib_deriv4_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_cdf.md)
body the Student t, the pseudo-Huber and the skew t register: the
components over the location and the scale from
[`loc_scale_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv_k.md),
and every component naming a shape parameter from
[`numerical_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv_k.md).

## Usage

``` r
partial_loc_scale_deriv_cdf_k(order)
```

## Arguments

- order:

  The derivative order, 3 or 4.

## Value

A function of `(distrib, q, theta, lower.tail, log, ...)` suitable for
registering as an S7 method on either generic, returning a named list of
numeric vectors of that order.

## Why the stencil is taken over everything

The stencil runs over every component and the closed ones then overwrite
it, so the location and the scale are computed twice. That is
deliberate: the whole cdf surface costs milliseconds at a thousand
quantiles, and the alternative, widening the stencil's signature to take
a subset at these orders, would touch a shared function for no
measurable gain. The orders below do use `which`, where the same
components are asked for far more often.

## Who is not here

The skew normal was among these families and is not any more. Owen's T
has elementary partial derivatives in both of its arguments, so its
shape components close too and it has a route of its own in
`cdf_skewnormal_higher.R`.

## See also

[`loc_scale_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv_k.md)
for the closed components;
[`numerical_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv_k.md)
for the rest;
[`partial_loc_scale_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_hess_cdf.md)
for the second order.
