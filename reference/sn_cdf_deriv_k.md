# Assemble the Skew Normal's CDF Derivatives of a Given Order

Chains the standard-coordinate table of
[`sn_cdf_std_derivs()`](https://statmodels7.github.io/distributions7/reference/sn_cdf_std_derivs.md)
through \\z = (q-\mu)/\sigma\\, the shape passing straight through as
the second index of that table. Together the two functions give the
family closed cdf derivatives at all four orders in all three
parameters.

## Usage

``` r
sn_cdf_deriv_k(distrib, q, theta, order)
```

## Arguments

- distrib:

  A `SkewNormal1Distrib` object, whose `params` name and order the
  components.

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu`, `sigma` (positive) and `alpha` (any
  sign).

- order:

  The derivative order, 1 to 4.

## Value

A named list of numeric vectors, derivatives of \\F\\ itself on the
natural scale, keyed as
[`deriv_names(distrib@params, order)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md).

## Details

The map is the one the other location-scale families are chained
through: \\z\\ is linear in the location and a reciprocal in the scale,
so a block naming the location twice contributes an exact zero and only
eight partials of the map are non-zero up to order four. The shape does
not enter \\z\\ at all, which is why it can be carried as an index of
the inner table rather than through the chain.

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale, \\\alpha\\ the
shape, \\z = (q-\mu)/\sigma\\ and \\F\\ the distribution function.

## See also

[`sn_cdf_std_derivs()`](https://statmodels7.github.io/distributions7/reference/sn_cdf_std_derivs.md)
for the inner table;
[`distrib_grad_cdf.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.SkewNormal1Distrib.md),
the family page;
[`loc_scale_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv_k.md)
for the same chain without a shape.
