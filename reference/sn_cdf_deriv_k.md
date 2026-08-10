# Assemble the Skew Normal's CDF Derivatives of a Given Order

Chains the standard-coordinate table of
[`sn_cdf_std_derivs`](https://statmodels7.github.io/distributions7/reference/sn_cdf_std_derivs.md)
through \\z = (q-\mu)/\sigma\\, the shape passing straight through as
the second index of that table.

## Usage

``` r
sn_cdf_deriv_k(distrib, q, theta, order)
```

## Arguments

- distrib:

  A `SkewNormal1Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu`, `sigma` and `alpha`.

- order:

  The derivative order, 1 to 4.

## Value

A named list of derivative components of \\F\\.

## See also

[`sn_cdf_std_derivs`](https://statmodels7.github.io/distributions7/reference/sn_cdf_std_derivs.md)
