# Higher Log-CDF Derivatives by the Compiled Quadrature

Builds the
[`distrib_deriv3_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_cdf.md)
or
[`distrib_deriv4_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_cdf.md)
body of the families whose compiled distribution function integrates its
shape derivatives and which are not location-scale: the derivatives of
\\F\\ of every order up to `order` from
[`cdf_tables()`](https://statmodels7.github.io/distributions7/reference/cdf_tables.md),
whose third and fourth orders are the exact integrals of
[`compiled_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/compiled_cdf_deriv_k.md),
put on the requested tail and scale by
[`cdf_scale_k()`](https://statmodels7.github.io/distributions7/reference/cdf_scale_k.md).

## Usage

``` r
compiled_deriv_cdf_k(order)
```

## Arguments

- order:

  The derivative order, 3 or 4.

## Value

A function of `(distrib, q, theta, lower.tail, log, ...)` suitable for
registering as an S7 method on either generic.
