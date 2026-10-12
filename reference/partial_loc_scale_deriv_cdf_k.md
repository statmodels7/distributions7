# Higher Log-CDF Derivatives When Only Some Parameters Are Location-Scale

Builds the
[`distrib_deriv3_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_cdf.md)
or
[`distrib_deriv4_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_cdf.md)
body that the two Student t, the two pseudo-Huber and the skew t
families register: the components over the location and the scale from
[`loc_scale_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv_k.md),
and every component naming a shape parameter from
[`continuous_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/continuous_cdf_deriv_k.md).

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

## The tables of lower order

At orders one and two every component comes from the family's own
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
and
[`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md),
and at orders three and four from
[`continuous_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/continuous_cdf_deriv_k.md),
which is the compiled quadrature for these families. The closed location
and scale components then replace the corresponding entries, so those
components are computed twice.

## The skew normal

The skew normal is not registered here. Owen's T has elementary partial
derivatives in both of its arguments, so the shape components of its
distribution function are closed, and it has a route of its own in
`cdf_skewnormal_higher.R`.

## See also

[`loc_scale_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv_k.md)
for the closed components;
[`numerical_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv_k.md)
for the rest;
[`compiled_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/compiled_cdf.md)
for the second order.
