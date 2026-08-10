# Higher Log-CDF Derivatives When Only Some Parameters Are Location-Scale

The higher-order companion of
[`partial_loc_scale_hess_cdf`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_hess_cdf.md):
the components over the location and the scale from the location-scale
construction, and the components naming a shape parameter from the
stencil.

## Usage

``` r
partial_loc_scale_deriv_cdf_k(order)
```

## Arguments

- order:

  The derivative order, 3 or 4.

## Value

A function suitable for registering as an S7 method.

## See also

[`loc_scale_cdf_deriv_k`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv_k.md)
